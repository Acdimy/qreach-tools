# TransitionSystem Optimization Record

## 2026-06-28: inQueue 数据结构优化

### 目标

本轮优化聚焦 `transition_system.hpp` 中显式朴素版 `qts_naive::TransitionSystem` 的 worklist 数据结构，优先服务实际工作流中更常用的 `postCondition` / `computingFixedPointPost()` 路径。

原实现中，重新入队时通过在 `currPostLocs` / `currPreLocs` 中执行 `std::find(...)` 判断 location 是否已经在队列中。该判断是线性扫描，在大规模 transition system 上可能被 fixed-point 迭代放大。

### 修改内容

在 `TransitionSystem` 中新增队列状态数组：

- `inPostQueue`
- `inPreQueue`

主要改动：

1. `setAnnotation(...)` 中避免重复向 `currPostLocs` / `currPreLocs` 入队。
2. `postConditionInit()` / `preConditionInit()` 中初始化并维护对应的 in-queue 状态。
3. `postConditionOneStep(...)` / `preConditionOneStep(...)` 中用 O(1) 的 `inPostQueue[loc]` / `inPreQueue[loc]` 检查替代原先的 `std::find(...)`。
4. `postConditions()` / `preConditions()` 中 location 出队后将对应 in-queue 标记恢复为 `false`。

核心效果是将 worklist 去重从 O(queue size) 降为 O(1)，同时更清楚地区分：

- `visitedPost` / `visitedPre`：是否曾经进入 fixed-point 传播；
- `inPostQueue` / `inPreQueue`：当前是否已经在 worklist 中。

### 构建结果

已重新构建 C++ library 和 Python extension：

```bash
cd python_pkg
../.venv/bin/python -m invoke build-qreach
../.venv/bin/python -m invoke build-pybind11
```

结果：构建成功。

### 单元测试结果

已运行以下测试：

```bash
../.venv/bin/python workflow_tests/test_simulation_grover.py
../.venv/bin/python workflow_tests/test_simulation_qft.py
../.venv/bin/python workflow_tests/test_ts_structure.py
```

结果：全部通过。

- `test_simulation_grover.py`：PASSED
- `test_simulation_qft.py`：PASSED
- `test_ts_structure.py`：All TransitionSystem structural tests PASSED

### benchmark 对比

#### 修改前

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.96 seconds
Time taken for model checking: 1.94 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 75.56 seconds
Time taken for model checking: 0.02 seconds
```

#### 修改后

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.99 seconds
Time taken for model checking: 1.89 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 76.82 seconds
Time taken for model checking: 0.02 seconds
```

### 结论

本轮修改保持了语义正确性，指定测试全部通过。`inQueue` 优化清理了 worklist 状态语义，并消除了 fixed-point 传播中的队列线性查找。

不过，从当前两个 benchmark 的结果看，性能提升并不明显，差异基本处于运行波动范围内。说明这些例子的主要瓶颈大概率不在 worklist 去重，而更可能在：

1. transition system 构造阶段；
2. Qiskit parsing / branch expansion；
3. `QOperation` 构造、复制、`postImage`、`disjunction`；
4. relation / location 数量本身；
5. 后续需要的 adjacency edge cache 或 lazy construction。

因此，本轮优化更适合作为后续 profiling、edge adjacency cache、lazy construction 的基础性低风险改动。

## 2026-06-29: post adjacency edge cache 优化

### 目标

本轮优化继续聚焦 `transition_system.hpp` 中显式朴素版 `qts_naive::TransitionSystem` 的 `postCondition` / `computingFixedPointPost()` 路径。

原 `postConditionOneStep(...)` 在遍历 `Locations[loc].postLocations` 时，会在热路径中反复通过：

```cpp
this->relations[std::make_tuple(loc, postLoc->idx)]
```

查找 relation。`relations` 当前是 `std::map<std::tuple<unsigned int, unsigned int>, QOperation>`，每次访问都需要构造 tuple 并做 map 查找；同时 `operator[]` 还有在 key 不存在时插入默认值的语义，不适合作为 fixed-point 传播热路径的只读访问方式。

### 修改内容

在 `TransitionSystem` 中新增 post adjacency edge cache：

```cpp
struct PostEdge {
    unsigned int to;
    const QOperation* relation;
};

std::vector<std::vector<PostEdge>> postEdgeCache;
void rebuildPostEdgeCache();
```

主要改动：

1. `addRelation(...)` 在新增 relation 后清空 `postEdgeCache`，避免缓存过期。
2. 新增 `rebuildPostEdgeCache()`，按 location 构造连续的后继边列表，每条边缓存：
   - successor location id；
   - 指向 `relations` 中对应 `QOperation` 的指针。
3. `postConditionInit()` 在 fixed-point post 计算开始前，如果缓存不存在或大小不匹配，则重建缓存。
4. `postConditionOneStep(...)` 改为遍历 `postEdgeCache[loc]`，并通过缓存中的 `const QOperation& relation` 访问 relation。
5. 对 post 热路径消除了重复的 `std::map` relation 查找和 repeated tuple-key `operator[]` 访问。

本轮只优化 post adjacency 路径，没有改变 pre-condition 路径，也没有改变 transition-system 构造策略；默认仍然是 eager construction。

### 构建结果

已重新构建 C++ library 和 Python extension：

```bash
cd python_pkg
../.venv/bin/python -m invoke build-qreach
../.venv/bin/python -m invoke build-pybind11
```

结果：构建成功。

备注：IDE/clang diagnostics 仍然报告本地 clangd include/config 相关问题，例如 Boost include 缺失和 STL 模板诊断；按照项目说明，以实际 `invoke build-qreach` / `invoke build-pybind11` 构建结果为准。

### 单元测试结果

已运行以下测试：

```bash
../.venv/bin/python workflow_tests/test_simulation_grover.py
../.venv/bin/python workflow_tests/test_simulation_qft.py
../.venv/bin/python workflow_tests/test_ts_structure.py
```

结果：全部通过。

- `test_simulation_grover.py`：PASSED
- `test_simulation_qft.py`：PASSED
- `test_ts_structure.py`：All TransitionSystem structural tests PASSED

### benchmark 对比

本轮 baseline 使用上一轮 `inQueue` 优化后的当前实现，在加入 post adjacency edge cache 前重新运行得到。

#### 修改前（inQueue-only baseline）

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.97 seconds
Time taken for model checking: 1.90 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 76.11 seconds
Time taken for model checking: 0.02 seconds
```

#### 修改后（post adjacency edge cache）

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.97 seconds
Time taken for model checking: 1.90 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 75.27 seconds
Time taken for model checking: 0.02 seconds
```

### 结论

本轮修改保持了语义正确性，指定测试全部通过。

post adjacency edge cache 清理了 `postConditionOneStep(...)` 的 relation 访问方式，使 post fixed-point 热路径不再依赖 `relations[tuple]` 的重复 map lookup，也避免了 `operator[]` 的隐式插入语义。这是后续继续优化 post 路径时更合适的数据结构基础。

从 benchmark 看：

1. `vqss` 的 build/check 时间基本不变；
2. `bv_n14` 的 build time 从 76.11s 到 75.27s，有小幅改善，但仍应视为运行波动范围内的弱信号；
3. 两个 benchmark 的 `model checking` / post-propagation 打印时间都很短或基本不变，说明当前主要瓶颈仍然不在这一次优化覆盖的 relation lookup 热点上。

后续更值得优先排查的方向仍然是：

1. transition system 构造阶段和连续测量导致的 location/relation 数量膨胀；
2. Qiskit parsing / branch expansion；
3. `QOperation` 构造、复制、`postImage`、`disjunction`；
4. 增加更细粒度的 profiling 以区分构造、fixed-point post、label/model-checking 各阶段耗时；
5. 进一步设计 lazy construction，尤其是从初始 annotation 前向传播的 workflow。
