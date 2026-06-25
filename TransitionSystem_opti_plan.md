# TransitionSystem Optimization Plan

## 背景与目标

当前 QReach 的近期重点是优化显式朴素版本 `qts_naive::TransitionSystem`（定义于 `transition_system.hpp`），而不是 QADD / symbolic transition system 路线。

显式 `TransitionSystem` 是当前 Qiskit 程序调试、annotation、labelling、qCTL / CTL workflow 的语义基准。实际使用中，`postCondition` / `computingFixedPointPost()` 往往比 `preCondition` 更重要，也更常被单独调用。因此，优化工作应优先围绕 forward post-condition fixed-point 展开，在保持现有 Python API 和语义兼容的前提下，降低时间和空间开销。

主要目标：

1. 优先优化 `postCondition` 热路径。
2. 降低大规模显式 transition system 中 location、relation、worklist 操作的开销。
3. 面向连续测量导致的指数级状态爆炸，设计后续 lazy construction 机制。
4. 使用现有 workflow 测试脚本评估优化效果和语义正确性。

## 当前朴素实现概览

`transition_system.hpp` 中的显式 transition system 主要由两个类构成：

- `qts_naive::Location`
- `qts_naive::TransitionSystem`

每个 `Location` 保存：

- 位置编号 `idx`
- 量子比特数 `qNum`
- `upperBound`
- `lowerBound`
- 前驱位置 `preLocations`
- 后继位置 `postLocations`
- classical proposition `cp`
- AP label 列表 `APs`

每条边的 quantum operation 保存在：

```cpp
std::map<std::tuple<unsigned int, unsigned int>, QOperation> relations;
```

其中 key 是 `(from, to)`。

### Post-condition fixed point

`computingFixedPointPost()` 的核心流程是：

1. `postConditionInit()` 初始化 worklist 和缓存。
2. `postConditions()` 循环弹出 `currPostLocs`。
3. `postConditionOneStep(loc)` 遍历 `loc` 的所有 successor。
4. 对每条边计算：

```cpp
QOperation postImage = Locations[loc].lowerBound.postImage(relations[(loc, postLoc)]);
postLoc->lowerBound = postLoc->lowerBound.disjunction(postImage);
```

5. 如果 successor 的 `lowerBound` 维度增长，则重新入队。

该流程是当前可达性分析和 Qiskit workflow 中最重要的路径。

## 已提出的优化点

### 1. 优先优化 postCondition 模块

实际 workflow 通常是：

```python
set_initial_state(ts, "000...")
ts.computingFixedPointPost()
```

也就是从程序开头 annotation 初始量子态，然后向前传播可达子空间。`preCondition` 在当前工作流中使用频率较低，因此第一阶段优化应集中在：

- `postConditionInit()`
- `postConditions()`
- `postConditionOneStep()`
- 与 post fixed-point 相关的 relation lookup、worklist、cache、dimension update 判断

`preCondition` 可以暂时保持语义稳定，只在结构改动自然影响时做最小兼容调整。

### 2. 使用 in-queue 标记替代 worklist 线性查找

当前代码在重新入队时使用：

```cpp
std::find(currPostLocs.begin(), currPostLocs.end(), postLoc->idx)
```

这会导致 worklist 中的线性扫描。对于 1600、16000 级别 location 的 transition system，这类扫描会随着 fixed-point 迭代放大。

建议增加明确的队列状态，例如：

```cpp
std::vector<bool> inPostQueue;
std::vector<bool> inPreQueue;
```

语义：

- 入队时设置 `inPostQueue[loc] = true`
- 出队时设置 `inPostQueue[loc] = false`
- bound 变化后，如果不在队列中，则重新入队

这样可以将队列去重从 O(queue size) 降为 O(1)。

### 3. 区分 visited 和 inQueue

当前 `visitedPost` 表示曾经访问过，但并不表示当前是否在队列里。因此代码仍需要扫描 `currPostLocs` 来避免重复入队。

建议将状态拆分为：

- `visitedPost` 或 `seenPost`：该 location 是否曾经被加入过 fixed-point 传播。
- `inPostQueue`：该 location 当前是否在 worklist 中。

这样逻辑更清晰，也更利于后续 profiling 和 lazy construction。

### 4. 优化 relation lookup

当前 relation 存储为：

```cpp
std::map<std::tuple<unsigned int, unsigned int>, QOperation> relations;
```

在 hot path 中访问方式为：

```cpp
relations[std::make_tuple(loc, postLocIdx)]
```

问题包括：

1. 每次构造 tuple。
2. `std::map` 是树结构查找，复杂度为 O(log E)。
3. `operator[]` 在 key 不存在时可能意外插入默认值。
4. transition traversal 已经通过 `postLocations` 找到 successor，却还需要再次查 map。

可选优化方向：

#### 方案 A：改为 unordered_map

将 `relations` 改为 hash map：

```cpp
std::unordered_map<std::tuple<unsigned int, unsigned int>, QOperation> relations;
```

优点：改动较小，保留现有 API 和数据结构。

缺点：仍然需要 tuple key lookup，仍然与 adjacency 分离。

#### 方案 B：在 adjacency 中直接保存 edge operation

增加显式 edge 结构：

```cpp
struct Edge {
    unsigned int from;
    unsigned int to;
    QOperation op;
};
```

或者分别保存：

```cpp
std::vector<std::vector<Edge>> outgoingEdges;
std::vector<std::vector<Edge>> incomingEdges;
```

这样 `postConditionOneStep()` 可以直接遍历：

```cpp
for (const auto& edge : outgoingEdges[loc]) {
    auto postLocIdx = edge.to;
    const QOperation& relation = edge.op;
    ...
}
```

优点：

- 避免 hot path 中反复 map lookup。
- 更适合 postCondition 优化。
- 后续 lazy construction 可以以 edge 为单位延迟生成。

缺点：

- 改动较大。
- 需要维护与现有 `relations`、`preLocations`、`postLocations` 的兼容。
- Python wrapper 或其他代码如果依赖 `getRelationName(from, to)`，仍需保留查询能力。

建议阶段性做法：

1. 第一阶段保留 `relations`，增加 adjacency edge cache。
2. `addRelation()` 同时维护旧结构和新结构。
3. post fixed-point 优先使用 `outgoingEdges`。
4. 对外 API 继续使用 `relations` 保持兼容。

### 5. 避免 identity self-loop 的多余计算

当前已有逻辑：

```cpp
if (postLoc->idx == loc && relations[(loc, postLoc->idx)].isIdentity) {
    continue;
}
```

如果引入 edge adjacency，可以把该判断移动到 edge 遍历中，并避免额外 relation lookup：

```cpp
if (edge.to == loc && edge.op.isIdentity) {
    continue;
}
```

也可以在构造阶段标记 identity self-loop，减少 fixed-point 热路径中的判断成本。

### 6. 审慎处理 computedTablePost

当前 post cache key 是：

```cpp
(loc, source_lower_dim, postLoc, target_lower_dim)
```

它只根据维度判断是否已经计算，而不是比较完整 `QOperation` 内容。

优点：

- 快。
- key 小。
- 当前语义和测试已经基于这个行为。

潜在风险：

- 不同 quantum subspace 可能维度相同但内容不同。
- 如果 lowerBound 内容变化但维度不变，缓存可能跳过一次实际需要的传播。

因此短期优化不建议随意改变该语义。可以先做 profiling 和 instrumentation，观察是否存在“维度不变但内容变化”的实际情况。若后续要增强正确性，需要单独设计更强的 bound version / generation counter，而不是简单替换。

### 7. 增加 profiling 统计

为了区分瓶颈来自 C++ 容器开销还是 `QOperation` 代数操作，应增加轻量 profiling。建议统计：

- post fixed-point 总耗时。
- `postConditionOneStep()` 调用次数。
- 遍历 edge 数量。
- 实际执行 `QOperation::postImage()` 次数。
- 实际执行 `disjunction()` 次数。
- 因 computed table 跳过的次数。
- 因 identity self-loop 跳过的次数。
- 重新入队次数。
- 最大 queue size。
- location 数和 relation 数。

可以沿用项目中已有的环境变量风格，例如：

```bash
TS_PROFILE=1 ...
```

第一阶段 profiling 的目的不是改变语义，而是给后续优化提供基线。

## 用户补充的重要优化方向

### 1. postCondition 比 preCondition 更重要

后续实现优先级应明确调整为：

1. `computingFixedPointPost()` 性能和空间优化。
2. post-only workflow 的 build / run 体验。
3. preCondition 只做必要兼容和回归保护。

也就是说，如果某些结构调整同时影响 pre/post，可以先保证 post path 获得实质收益；pre path 不应成为第一阶段优化的主要复杂度来源。

### 2. 连续测量导致显式构造状态爆炸

量子程序中如果存在多个 qubit 的连续 measurement，即使实际 forward reachability 从给定初始态出发最终只会走少数分支，当前 parser / transition system construction 仍然会预先构造所有分支位置。

对于连续测量，分支数量可能指数增长。例如测量 k 个 qubit 时，显式控制流可能产生接近 `2^k` 个分支位置。每个状态之间还会生成对应 relation quantum operation，导致：

- location 数指数增长；
- relation 数指数增长；
- relation 上的 `QOperation` 占用大量空间；
- 后续 post fixed-point 还要遍历大量实际不可达或无需展开的边。

这对大程序尤其明显。

### 3. Lazy construction 开关

后续可以为 explicit transition system 增加 lazy construction 模式。

适用场景：

- annotation 只在程序开头设置初始量子态；
- 不需要在未展开位置上提前设置 quantum annotation；
- 主要运行 forward `postCondition`；
- 连续 measurement 后的许多分支在当前 lowerBound 下不会实际到达。

理想行为：

1. parser 不一次性展开所有 measurement 分支后的 locations / relations。
2. 遇到连续测量区域时，先保存可展开的 construction thunk / continuation / deferred block。
3. `computingFixedPointPost()` 推导到该位置时，根据当前 location 的 `lowerBound` 和 measurement relation 判断实际可达分支。
4. 只为实际需要传播的分支生成后续 locations 和 relations。
5. 新生成 location 后继续参与 post worklist。

潜在收益：

- 显著降低连续 measurement 程序的构造期空间开销。
- 减少不可达分支 relation 的 `QOperation` 创建。
- 对 post-only workflow，减少无用 location 遍历。

需要注意：

- lazy construction 会影响 parser、metadata、marker、instruction location mapping。
- 如果用户需要完整 transition system 或需要对所有 location 做 annotation / labelling，则 lazy 模式可能不适用。
- 需要明确提供开关，默认可以先保持 eager construction，避免破坏现有语义。

可能 API：

```python
parse_qiskit_cir(..., lazy=False)
```

或 TransitionSystem 构造开关：

```python
ts = pyqreach.TransitionSystem(lazy=True)
```

具体放在 parser 还是 TS 层，需要进一步设计。由于 lazy construction 涉及 Qiskit control-flow lowering，最终可能需要 parser 和 TS 同时配合。

## 分阶段实施计划

### 阶段 1：建立 profiling 和测试基线

目标：不改变语义，先得到可比较的数据。

任务：

1. 在 post fixed-point 中增加可选 profiling。
2. 统计 location 数、relation 数、edge traversal、postImage、disjunction、skip、入队次数等。
3. 运行两个 workflow 测试脚本记录 baseline。
4. 保留当前输出格式，避免破坏现有脚本。

建议命令：

```bash
cd python_pkg
TS_PROFILE=1 ../.venv/bin/python workflow_tests/test_vqss_correct.py
TS_PROFILE=1 ../.venv/bin/python workflow_tests/test_bv_n14.py
```

### 阶段 2：worklist O(1) 去重

目标：低风险优化 postCondition worklist。

任务：

1. 增加 `inPostQueue`。
2. 在 post 入队 / 出队时维护该状态。
3. 用 O(1) 检查替代 `std::find(currPostLocs.begin(), currPostLocs.end(), ...)`。
4. 对 preCondition 可选择同步增加 `inPreQueue`，但不是优先目标。
5. 跑测试脚本确认结果不变。

风险较低，适合作为第一批实际代码改动。

### 阶段 3：post adjacency edge cache

目标：减少 post fixed-point 热路径中的 relation map lookup。

任务：

1. 增加 edge 结构或 outgoing adjacency cache。
2. `addRelation()` 同时维护：
   - `relations[(from, to)]`
   - `Locations[from].postLocations`
   - `Locations[to].preLocations`
   - `outgoingEdges[from]`
   - 可选 `incomingEdges[to]`
3. `postConditionOneStep()` 改为遍历 `outgoingEdges[loc]`。
4. 保留 `relations` 用于兼容现有 API，例如 `getRelationName(from, to)`。
5. 跑测试脚本对比 location 数、结果、耗时。

### 阶段 4：relation container 兼容性优化

目标：进一步减少 relation 查询和内存开销。

可选任务：

1. 将 `relations` 从 `std::map` 改为 `std::unordered_map`。
2. 避免 hot path 使用 `operator[]`。
3. 对只读查询使用 `find()` / `at()`。
4. 如果 edge adjacency 已覆盖 hot path，则该阶段收益可能较小，可根据 profiling 决定是否实施。

### 阶段 5：lazy construction 设计原型

目标：解决连续 measurement 造成的空间爆炸。

任务：

1. 分析 Qiskit parser 中 measurement / if_test / control-flow lowering 的 construction 点。
2. 定义 lazy 模式下哪些结构可以延迟：
   - location
   - relation
   - branch continuation
   - marker metadata
3. 定义何时 materialize lazy branch：
   - post propagation 到达 measurement source location；
   - 当前 lowerBound 对某 measurement outcome 有非零 postImage；
   - 需要访问某 marker / instruction location。
4. 明确 lazy 模式的限制：
   - 不适合需要完整 TS 的工作流；
   - 不适合在任意未展开位置预设 annotation；
   - metadata 可能需要 lazy-aware representation。
5. 先选一个连续 measurement benchmark 做原型，不直接大范围重构 parser。

## 推荐测试脚本

后续优化以以下两个脚本作为主要 workflow benchmark。

### 1. `python_pkg/workflow_tests/test_vqss_correct.py`

用途：Verifiable quantum secret sharing workflow。

特点：

- 构造 14-qubit QReachCircuit。
- 包含 Steane code preparation。
- 使用 inline marker：`enc`。
- 包含多个 measurement 和 classical `if_test`。
- 执行完整 workflow：
  - parse Qiskit circuit；
  - set initial state；
  - `computingFixedPointPost()`；
  - `label_snapshot()`；
  - `annotate(ts, ["leaf"])`；
  - `modelChecking(ts, 'AG (leaf -> target)')`。

运行命令：

```bash
cd python_pkg
../.venv/bin/python workflow_tests/test_vqss_correct.py
```

profiling 命令：

```bash
cd python_pkg
TS_PROFILE=1 ../.venv/bin/python workflow_tests/test_vqss_correct.py
```

当前参考规模：

- 运行时间约 1 秒级。
- location 数约 1600 级别。

优化后应检查：

- `Transition System Locations` 是否保持一致，除非启用 lazy 模式。
- model checking result 是否保持一致。
- build time 和 post/model checking time 是否改善。

### 2. `python_pkg/workflow_tests/test_bv_n14.py`

用途：Bernstein-Vazirani `bv_n14.qasm` benchmark workflow。

特点：

- 从 QASM 文件读取 circuit：

```python
benchmark/benchpress-medium/supported/bv_n14.qasm
```

- parse Qiskit circuit；
- set all-zero initial state；
- 执行 `computingFixedPointPost()`；
- 打印 location 数、build time、model checking time。

运行命令：

```bash
cd python_pkg
../.venv/bin/python workflow_tests/test_bv_n14.py
```

profiling 命令：

```bash
cd python_pkg
TS_PROFILE=1 ../.venv/bin/python workflow_tests/test_bv_n14.py
```

当前参考规模：

- 运行时间约 60 秒级。
- location 数约 16000 级别。

该脚本更适合作为性能优化主 benchmark。优化后应重点比较：

- build time；
- `computingFixedPointPost()` time；
- location 数；
- postImage 次数；
- disjunction 次数；
- relation lookup / edge traversal 次数；
- peak queue size。

## 正确性与兼容性要求

优化过程中应保持以下行为：

1. 现有 Python API 不应破坏：
   - `pyqreach.TransitionSystem()`
   - `parse_qiskit_cir(...)`
   - `set_initial_state(...)`
   - `computingFixedPointPost()`
   - `label_snapshot(...)`
   - `annotate(...)`
2. 默认模式下仍应 eager 构造完整 transition system。
3. `getLocationNum()` 在非 lazy 模式下应保持一致。
4. `getRelationName(from, to)` 等 relation 查询接口应保持可用。
5. `reached` / `valid` annotation 语义仍应基于：

```python
ts.printDims(loc)[1] > 0
```

6. 不应将当前优化转向 QADD / symbolic transition system。
7. 不应在没有明确设计的情况下改变 `computedTablePost` 的语义。

## 初步优先级

建议优先级如下：

1. Profiling baseline。
2. Post worklist `inPostQueue` 优化。
3. Post adjacency edge cache，减少 relation map lookup。
4. Relation container 清理和 `operator[]` 避免。
5. Lazy construction 设计与原型。
6. 根据 profiling 决定是否进一步深入 `QOperation::postImage` / `disjunction` 层。

## 备注

当前优化方向应围绕显式朴素 `TransitionSystem`，尤其是 Qiskit workflow 中的 forward reachability。QADD / symbolic implementation 暂不作为本计划的一部分。
