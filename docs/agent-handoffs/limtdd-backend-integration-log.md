# LimTDD 后端替换 — QReach 侧集成日志与已知问题

> **Audience:** 接续后端替换工作的 agent（QReach 侧 + LimTDD 侧）。
> **Status:** 集成进行中 —— Phase 0–5 已完成、CFLOBDD 封装完毕、LimTDD 已接线；剩余 dense 实现与若干 qNum/level 约定待办。
> **Last updated:** 2026-08-28（覆盖 08-13 → 08-17 的进度与发现，此前文档未收录）。
> **分支:** `limtdd-backend`。
> **姊妹文档:** `limtdd-backend-implementation-status.md`（LimTDD 侧实现现状）、`limtdd-backend-qreach-plan.md`（QReach 侧计划）、`backend-replacement-api-contract.md`（API 契约）。

本文档补记 `limtdd-backend-qreach-plan.md`（停在 Phase 计划）与 `limtdd-backend-implementation-status.md`（停在 08-13 LimTDD 侧现状）之后、尚未落入仓库文档的进度与发现。

---

## 1. 集成进度：Phase 0–5 全部完成

`docs/agent-handoffs/limtdd-backend-qreach-plan.md` 定义的 Phase 0–6 中，Phase 0–5 已全部落地并验证，Phase 6（full verification + 正式交接）部分完成。各阶段结果：

| Phase | 内容 | 结果 |
|---|---|---|
| 0 | `dd_backend.hpp` 接口头 + `cflobdd/CFLOBDD/dd_backend_cflobdd.h` CFLOBDD 透传实现 | 完成。CFLOBDD 初始化是 **7 次调用**：3× `CFLOBDDNodeHandle::Init*` + `InitPairProductCache()` + `InitTripleProductCache()` + `Matrix1234Initializer()` + `VectorInitializer()`；漏掉 product caches 会在首次 `KroneckerProduct`/`MatrixMultiply` 时于 `Hashtable::Fetch` 段错误 |
| 1 | `quantum_operation.hpp` 换 include/类型别名（`CFLOBDD_COMPLEX_BIG`→`DD`、`BIG_COMPLEX_FLOAT`→`DDComplex`） | 完成 |
| 2 | 直接 DAG 访问替换为查询函数（`checkifzero`→`IsApproximatelyZero`、`.root->level`→`GetLevel`、`dot`→`InnerProduct`、`normalize`→`Normalize`）；transpose 损坏 fallback 下沉到 backend | 完成 |
| 3 | 删除自由函数 `ApplyGateF*`/`InitializeWithVector`，改调 `DDMatrix::MkSingleQubitGateOnN*`/`DDVector::InitializeWithAmplitudes` | 完成 |
| 4 | `concretize()` 全部 `Matrix1234ComplexFloatBoost::Mk*`→`DDMatrix::Mk*`、`VectorComplexFloatBoost::`→`DDVector::`；`MatrixMultiplyV4WithInfo` 拆分为 `MatrixMultiply`(matrix×matrix) 与 `MatrixMultiplyWithVector`(gate×vector) | 完成。`quantum_operation.hpp` 达到 **100% backend-clean**（grep 无 `CFLOBDD`/`BIG_COMPLEX`/`CFL_OBDD`/`.root->`/`SH_OBDD`） |
| 5 | `transition_system.hpp` 的 7 调用初始化替换为 `DDMatrix::Initialize(); DDVector::Initialize();`（先 matrix 后 vector） | 完成。CFLOBDD 已完整封装进 `dd_backend.hpp` + `dd_backend_cflobdd.h` |

**LimTDD 接线（`QREACH_USE_LIMTDD`）**：
- `dd_backend.hpp` 的 LimTDD 分支 include `dd/backend/DDVector.hpp` + `dd/backend/DDMatrix.hpp`（LimTDD 在**全局命名空间**提供完整内联实现），别名 `qreach::DD = limtdd::DD`、`qreach::DDComplex = limtdd::DDComplex`。DDVector/DDMatrix 的前向声明只在 CFLOBDD 分支里。
- Makefile：`USE_LIMTDD=1` → `-DQREACH_USE_LIMTDD` + `-I$(LIMTDD_PATH) -I$(LIMTDD_INCLUDE)`，并编译 3 个核心 `.cpp`（`dd/{Edge,Maps,Node}.cpp`）**取代** CFLOBDD 源文件。`LIMTDD_PATH=../LimTDDexpr/LimTDD/DDPackage`、`LIMTDD_INCLUDE=../LimTDDexpr/include`（xtensor/xtl）。
- 冒烟：`scripts/dd_backend_smoke.cpp`（backend primitives）与 `scripts/limtdd_explicit_smoke.cpp`（Bell H+CX 后像，dim [4,1]）均通过。
- ⚠️ `USE_LIMTDD=1 make all` 会把 LimTDD 的 `.o` 写进外部项目目录 `../LimTDDexpr/.../dd/`（污染外部仓库）；切回 CFLOBDD 需 `rm libqreach.so && make all`。

---

## 2. 已解决的后端对接 bug（08-13 → 08-17）

按时间顺序，均为跨后端对接暴露、现已修复/缓解：

### 2.1 `printFormal` 后端格式耦合（08-13，已修）

`quantum_operation.hpp::printFormal()` 原先解析 `VectorPrintColumnHead` 的字符串输出（tokenize + `std::stod`/`sscanf`），隐含 CFLOBDD 的 dense `(re,im)` 格式。LimTDD 的 `VectorPrintColumnHead` 输出 `|idx>\tre im` 行，导致 `stod` 抛 "no conversion"。

**修复**：`printFormal` 改为直接枚举 `DDVector::GetNonZeroAmplitudes`（backend-neutral）。CFLOBDD 输出不变（`test_RUS`/`test_newapi` 在 CFLOBDD 下 exit 0 验证）。

### 2.2 `InnerProduct` scaling bug（08-14，已修，LimTDD 侧）

详见 `limtdd-innerproduct-bug-report.md`。`DDVector::InnerProduct` 对 ≥2-qubit 向量返回错误标量：`<basis|basis>` 被缩放 `2^(n-1)`（1q→1、2q→2、4q→8），叠加项还有**非均匀**因子（`<v1|v3>` = -24× = -3·2³）。因缩放非均匀，Gram-Schmidt 的 `⟨vi|v⟩/⟨vi|vi⟩` 比值被破坏，线性相关向量无法坍缩 → RUS workflow 多出支撑向量（CFLOBDD 2 vs LimTDD 4）。

**修复**：LimTDD 侧将 `InnerProduct` 改为基于 `GetNonZeroAmplitudes(a,0)` × `GetNonZeroAmplitudes(b,0)`（按大端索引 join）的 sparse dot product，绕过 `cont` 的坏标量路径。验证：最小复现 8 项全过（`<basis|basis>`=1 在 1/2/4q；投影系数精确）；`span_qops([v1..v4])`=2；RUS L5=3、L15–L19=2，与 CFLOBDD 一致；`test_RUS`/`newapi`/`grover`/`lazy_measurement` 在 LimTDD 下全绿。

### 2.3 `GetLevel` 约定不匹配（08-14 识别 → 08-17 修）

`SingleVecTerm(DD x)`（`quantum_operation.hpp`）原先设 `qNum = 2^(GetLevel(x)-1)`，隐含 CFLOBDD 的 matrix-form 存储（`VectorToMatrixInterleaved` 把向量 level 抬到 matrix level，故 `2^(level-1)`=qubits）。LimTDD 的 `VectorToMatrixInterleaved` 是 no-op 且 `GetLevel` 返回**向量** level，导致 qNum 减半。症状：`examples/test_rus_wp.py` 在 `projectIn` 崩溃（`oplist[0]->qNum == vec.qNum` 断言失败，4 vs 2）；`grover_wp` 结果支撑向量数 8 vs CFLOBDD 3。

**修复**：引入 backend-neutral 的 `QubitsFromLevel(level)`（见 §2.4），`SingleVecTerm(DD)` 改为 `qNum = QubitsFromLevel(GetLevel(x))`。

### 2.4 Backend-neutral padding/level 重构（08-17，已落地）

power-of-2 padding（`qNum = 2^ceil(log2(n))`）是 CFLOBDD 特有，LimTDD 原生表示任意 qubit 数。在 `dd_backend.hpp` 的 `qreach` 命名空间（两个分支都定义）新增：

```cpp
NormalizeQubitCount(q)     // CFLOBDD: 向上取 2 的幂；LimTDD: q
VectorLevelForQubits(q)    // CFLOBDD: log2(q)；LimTDD: q
MatrixLevelForQubits(q)    // CFLOBDD: log2(q)+1；LimTDD: q
QubitsFromLevel(level)     // CFLOBDD: 2^(level-1)；LimTDD: 2^level
```

`quantum_operation.hpp` 的所有 level/qNum 运算改走这些 helper（gate ctor、`SingleVecTerm` string/amp/DD 三个 ctor、`concretize` 的 state_level/n、`projectIn`、`resetall`、`genProjMeasSpace`、自由函数 `logicqNum`），并移除了 power-of-2 断言。

- **CFLOBDD 等价性已验证**：VQSS correct→True / buggy→False，`test_ts_structure`/`test_simulation_grover` 通过。
- **LimTDD 侧需配套改**（见 §3.2），否则语义仍不完整。

---

## 3. 仍待办 / 已识别的剩余问题

### 3.1 LimTDD 侧待办

1. **dense O(4^n) 实现**（已在 `limtdd-backend-implementation-status.md` §4 记录）：`MkSingleQubitGateOnN`/`MkCNOT` 等多 qubit 门、`KroneckerProduct`、`MatrixMultiply`、`Conjugate`、`Transpose` 都是「构建 dense 矩阵 → `to_tdd`」或「枚举 → 变换 → 重建」，只适用 n≤10–12。实测 `grover_wp` 比 CFLOBDD 慢约 **2300×**（23s vs 0.01s）。根因是 `cont` 的两个问题：变量序必须交错（已用交错预注册解决），以及对 disjoint 矩阵 Kron 后再 `cont` 产生畸形 DD（路径变量号重复 `v=0→v=0`），故 Kron/门提升不能走「Kron(cont) + 再 cont」。

2. **`kMaxDenseQubits=12` 与 QReach padding 冲突**：10/14-qubit VQSS 例子 pad 到 16 后，`denseSingleQubitGateOnN` 抛 "too many qubits"。CFLOBDD 能跑（correct→True、buggy→False），LimTDD 不能。与 §2.4 的 `NormalizeQubitCount` 语义相关，需 LimTDD 侧决定如何表示非 power-of-2 qubit 数。

3. **门构造器 level→n 约定待改**（§2.4 重构后新暴露）：LimTDD 侧以下函数当前按 `level`（log2 尺度）解释参数，需改为按 `n`（qubit 数）解释，才能配合 `VectorLevelForQubits(q)=q`（LimTDD 分支）：
   - `MkBasisVector`、`NoDistinctionNode`、`MkSwap`、`MkiSwap`、`MkCP`（当前内部用 `qubitCount(level)=2^level` / `n=2^(level-1)`）；
   - `MkCNOT`/`MkCCNOT`/`MkSingleQubitGateOnN` 已按 `n`，无需改。
   - `GetLevel` 当前返回 `ceil(log2(n))`，对非 power-of-2 是**有损**的；非 power-of-2 的 `SingleVecTerm(DD)` 需要 `GetLevel` 返回 `n`（或另加 `GetQubitCount`）。

### 3.2 QReach 侧待办

4. **`dd_backend_cflobdd.h` 的 header-scope `using namespace CFL_OBDD;`**：命名空间污染，应收敛进 `DDVector`/`DDMatrix` 命名空间内。

5. **`transition_system_qadd.hpp`（SHELVED 符号化路径）**：仍直接裸调 CFLOBDD init（通过传递 using 仍能工作），且该路径是 CFLOBDD-only，`test_qreach`/wrapper 在 `QREACH_USE_LIMTDD` 下无法 build —— 需加 guard（或明确其 CFLOBDD-only 身份）。

### 3.3 外观性差异（不影响正确性）

6. **`GetNonZeroAmplitudes` 索引约定不一致**：CFLOBDD 返回 returnMap **位置**（非二进制值，故 `printFormal` 的 bitstring 标签"错误"），LimTDD 返回**大端二进制值**。振幅/子空间一致，仅 bitstring 标签不同。

---

## 4. 验证清单（当前状态）

**静态（backend-clean 完成检查）**：
```bash
grep -nE "Matrix1234ComplexFloatBoost|VectorComplexFloatBoost|CFLOBDD|BIG_COMPLEX|CFL_OBDD|CFLOBDDNodeHandle|\.root->|SH_OBDD" \
  quantum_operation.hpp transition_system.hpp cl_proposition.hpp python_pkg/qreach_python_wrapper.cpp
# 期望：无匹配
```

**CFLOBDD（默认后端，每阶段 gate）**：
```bash
make clean && make test && ./test_qreach 8
cd python_pkg && ../.venv/bin/python -m invoke build-qreach && ../.venv/bin/python -m invoke build-pybind11
../.venv/bin/python workflow_tests/test_newapi.py
../.venv/bin/python workflow_tests/test_grover.py
../.venv/bin/python workflow_tests/test_RUS.py
../.venv/bin/python workflow_tests/test_lazy_measurement.py
../.venv/bin/python workflow_tests/symbolic/test_symts_minimal.py   # 符号化回归（shelved，仍须过）
```

**LimTDD（冒烟）**：
```bash
USE_LIMTDD=1 BOOST_PATH=../BOOST/boost_1_81_0 make test
# scripts/dd_backend_smoke.cpp / scripts/limtdd_explicit_smoke.cpp
```

---

## 5. 相关文档

- `backend-replacement-api-contract.md` — Document A（API 契约，LimTDD 侧实现依据）。
- `limtdd-backend-qreach-plan.md` — QReach 侧计划（Phase 0–6 定义，本文档记录其执行结果）。
- `limtdd-backend-implementation-status.md` — LimTDD 侧实现现状（08-13）。
- `limtdd-backend-implementation-brief.md` — 交给 LimTDD agent 的一页简报。
- `limtdd-innerproduct-bug-report.md` — InnerProduct scaling bug 详情与 resolution。
- `cflobdd-transpose-dag-corruption.md`、`cflobdd-level8-transpose-bug.md` — CFLOBDD transpose/level≥7 bug 背景（dot/normalize fallback 的由来）。
