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
