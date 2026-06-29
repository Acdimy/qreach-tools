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

## 2026-06-29: relation container unordered_map 兼容性优化

### 目标

本轮优化针对 `transition_system.hpp` 中显式朴素版 `qts_naive::TransitionSystem` 的 relation 容器本身。

此前 `relations` 使用：

```cpp
std::map<std::tuple<unsigned int, unsigned int>, QOperation> relations;
```

这会让 relation 查询和插入走 ordered map 的树结构。上一轮 `post adjacency edge cache` 已经减少了 post fixed-point 热路径中的 repeated lookup；本轮进一步把底层 relation container 切换为 hash map，并清理只读查询中的 `operator[]` 用法，避免只读路径出现隐式插入语义。

### 修改内容

1. 为二元 relation key 增加 hash specialization：

```cpp
std::hash<std::tuple<unsigned int, unsigned int>>
```

2. 将 `TransitionSystem::relations` 从 `std::map` 改为：

```cpp
std::unordered_map<std::tuple<unsigned int, unsigned int>, QOperation> relations;
```

3. 保持 `relations` 字段名和 pybind 暴露方式不变，因此 Python 侧现有访问方式（例如 `len(ts.relations)`）保持兼容。

4. `getRelationName(...)` 和 `rebuildPostEdgeCache()` 继续使用 `find()` 做只读查询。

5. 清理 pre-condition 路径中的只读 relation 访问，将多处：

```cpp
this->relations[std::make_tuple(...)]
```

替换为一次 `find()` 后复用：

```cpp
auto relationIt = this->relations.find(std::make_tuple(...));
assert(relationIt != this->relations.end());
const QOperation& relation = relationIt->second;
```

6. `addRelation(...)` 仍然作为写入路径使用 `operator[]` 更新 relation，这是有意保留的写入语义；只读 fixed-point 查询路径不再依赖 `operator[]`。

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

本轮 baseline 使用上一轮 `post adjacency edge cache` 优化后的当前实现，在改为 `std::unordered_map` 前重新运行得到。

#### 修改前（post adjacency edge cache baseline）

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.95 seconds
Time taken for model checking: 1.91 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 72.71 seconds
Time taken for model checking: 0.02 seconds
```

#### 修改后（relations unordered_map）

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.93 seconds
Time taken for model checking: 1.84 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 72.13 seconds
Time taken for model checking: 0.02 seconds
```

### 结论

本轮修改保持了语义正确性，指定测试全部通过。

relation container 改为 `std::unordered_map` 后，relation 的平均查找/更新复杂度从 ordered map 的 O(log E) 变为 hash map 的平均 O(1)。同时，只读 fixed-point 查询路径改用 `find()` 后复用 `const QOperation&`，避免了 `operator[]` 在只读路径上的隐式插入语义，也与 `TransitionSystem_opti_plan.md` 中“对只读查询使用 `find()` / `at()`”的方向一致。

从 benchmark 看：

1. `vqss` 从 0.95s build / 1.91s check 变为 0.93s build / 1.84s check，有小幅改善；
2. `bv_n14` build time 从 72.71s 到 72.13s，有小幅改善；
3. 改善幅度仍然不大，应视为低风险数据结构优化和兼容性清理，而不是主要瓶颈已解决。

后续如果继续优化，建议优先增加细粒度 profiling，确认时间主要消耗在：

1. parser / transition-system 构造；
2. relation/location 生成数量；
3. `QOperation` 复制和运算；
4. post fixed-point propagation；
5. labelling / model checking 输出阶段。

## 2026-06-29: lazy measurement construction 初步实现

### 目标

本轮优化针对连续 measurement 在 eager transition-system 构造阶段造成的指数级分支扩张。根据当前语义约束，本轮不把 transition system 改成“天然 reachable-only graph”，而是只实现 measurement-focused lazy construction：

1. measurement 的两个 outcome location 仍然显式构造；
2. quantum post-image 非零的 outcome 继续进入后续 parser frontier；
3. quantum post-image 为零的 outcome 构造为 leaf placeholder，保留 classical AP 和 incoming measurement edge；
4. zero-outcome placeholder 不再继续展开后续控制流；
5. 由用户显式调用 `parse_qiskit_cir_lazy(...)` 开启，默认 `parse_qiskit_cir(...)` eager 行为保持兼容。

这样保留了现有 CTL 使用方式中显式 `reached` / `valid` AP 的语义：不带 reach guard 的公式仍然能看到 graph 中的 zero-bound placeholder；带 reach guard 的公式则可以过滤这些 quantum-unreachable locations。

### 修改内容

1. 在 `pyqreach.QOperation` Python binding 中新增：

```python
op.post_image(relation)
op.dim()
op.is_zero()
```

用于 Python parser 在构造期间判断 measurement outcome 的 quantum reachability。

2. 扩展 `ParseResult` lazy metadata：

```python
lazy: bool
lazy_pruned_locations: list[int]
lazy_pruned_by_instruction: dict[int, list[int]]
```

3. 新增 `parse_qiskit_cir_lazy(...)` wrapper。该 wrapper 要求提供 `initial_state` 或 `initial_op`，因为 lazy measurement 需要在 parse 过程中前向维护 lowerBound。

4. 在 parser 内部增加 lazy context。lazy mode 下普通 gate / control-flow identity branch 会同步传播 lowerBound；measurement 分支则根据 post-image 区分：

- 非零 post-image：构造正常 successor，更新 lowerBound，加入下一轮 `currLoc`；
- 零 post-image：构造 placeholder successor，保留 classical AP 和 relation，记录到 `lazy_pruned_locations`，但不加入 `currLoc`。

5. `while_loop` 在 lazy mode 下暂时显式报 `NotImplementedError`，避免在尚未实现 lazy loop skeleton 的情况下产生不清晰语义。

6. 新增测试/benchmark 脚本：

```bash
python_pkg/workflow_tests/test_lazy_measurement.py
python_pkg/workflow_tests/test_vqss_correct_lazy.py
python_pkg/workflow_tests/test_bv_n14_lazy.py
```

### 构建结果

已重新构建 Python extension：

```bash
cd python_pkg
../.venv/bin/python -m invoke build-pybind11
```

结果：构建成功。

### 单元测试结果

已运行：

```bash
PYTHONPATH=. ../.venv/bin/python workflow_tests/test_lazy_measurement.py
PYTHONPATH=. ../.venv/bin/python workflow_tests/test_simulation_grover.py
PYTHONPATH=. ../.venv/bin/python workflow_tests/test_simulation_qft.py
PYTHONPATH=. ../.venv/bin/python workflow_tests/test_ts_structure.py
```

结果：全部通过。

- `test_lazy_measurement.py`：All lazy measurement tests PASSED
- `test_simulation_grover.py`：PASSED
- `test_simulation_qft.py`：PASSED
- `test_ts_structure.py`：All TransitionSystem structural tests PASSED

### benchmark 对比

#### eager baseline

`workflow_tests/test_vqss_correct.py`：

```text
Transition System Locations: 1635
Time taken for building transition system: 0.99 seconds
Time taken for model checking: 1.91 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14.py`：

```text
Transition System Locations: 16424
Time taken for building transition system: 75.84 seconds
Time taken for model checking: 0.02 seconds
```

#### lazy measurement

`workflow_tests/test_vqss_correct_lazy.py`：

```text
Transition System Locations: 315
Lazy pruned locations: 40
Time taken for building transition system: 0.20 seconds
Time taken for model checking: 0.30 seconds
Model checking result: True
```

`workflow_tests/test_bv_n14_lazy.py`：

```text
Transition System Locations: 68
Lazy pruned locations: 13
Time taken for building transition system: 0.03 seconds
Time taken for model checking: 0.02 seconds
```

### 结论

lazy measurement construction 对当前两个 benchmark 的构造阶段有明显改善：

1. `vqss` location 数量从 1635 降到 315，build time 从约 0.99s 降到约 0.20s；
2. `bv_n14` location 数量从 16424 降到 68，build time 从约 75.84s 降到约 0.03s；
3. zero-outcome measurement branch 仍然显式存在于 transition graph 中，并保留 classical AP / incoming edge；
4. pruning 只阻止 zero-bound placeholder 继续展开后续控制流，因此没有把现有 TS 语义改成天然 reachable-only。

本轮实现仍是第一阶段：重点支持 initial annotation / forward post workflow 下的 lazy measurement。后续如需完整兼容复杂 control flow 中 unreachable suffix 的 marker / annotation 精确位置，可以继续设计更细的 unreachable skeleton mode。