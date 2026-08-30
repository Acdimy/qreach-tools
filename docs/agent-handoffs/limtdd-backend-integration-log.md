# LimTDD 后端替换 — QReach 侧集成日志与已知问题

> **Audience:** 接续后端替换工作的 agent（QReach 侧 + LimTDD 侧）。
> **Status:** 集成**已重建并验证**（2026-08-30）——QReach 语义层 + 显式 TS 在 LimTDD 后端下全绿（含 VQSS）；LimTDD 侧已紧凑化门/态/Kron、修好 level→n 约定。剩余 `Conjugate`/`Transpose` dense 化与矩阵代数紧凑化。
> **Last updated:** 2026-08-30（覆盖 08-13 → 08-30；§5 记录 08-30 重同步结果）。
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

### 3.1 LimTDD 侧待办（2026-08-30 更新：多数已由 LimTDD 侧完成）

> ⚠️ 本节原列的三条（dense 门、`kMaxDenseQubits` 冲突、level→n 约定）**均已在 LimTDD `limtdd-backend` 分支（HEAD `77936fb`）上解决**。详见 §5 的重同步验证。以下为更新后的剩余项。

1. ~~**dense O(4^n) 门构造**~~ → **已紧凑化**：`MkSingleQubitGateOnN` 现走 `cont`（"no dense O(4^n)"），`KroneckerProduct` 也走紧凑张量积。`grover_wp` 的 2300× 性能差应已大幅改善（待重新测量确认）。

2. ~~**`kMaxDenseQubits=12` 与 padding 冲突**~~ → **已解除**：`kMaxDenseQubits=12` 现在只 guard dense 枚举辅助（`matrixToDense`），不再 guard 门构造。VQSS（10/14 qubit，pad 到 16）现已可跑（§6 验证：correct→True、buggy→violation）。

3. ~~**门构造器 level→n 约定**~~ → **已改**：`MkBasisVector`/`NoDistinctionNode`/`MkSwap`/`MkiSwap`/`MkCP` 现收 `n`（qubit 数）；`MkCNOT`/`MkCCNOT` 保留 `/*level*/`（未用）+ `n`；`GetLevel` 返回 `n`（= `tdd.index_set.size()`）。与 QReach 侧 `VectorLevelForQubits(q)=q` / `QubitsFromLevel(level)=level` 对齐。

4. **`Conjugate` / `Transpose` 仍是 dense 实现**（**剩余硬骨头**）：`DDMatrix.hpp` 中 `Conjugate` 标注 "Implemented via dense"（`matrixToDense` → 逐项共轭 → 重建）。这是 LimTDD 侧自认的下一步（HISTORY.md 指向 `limtdd-backend-conjugate-transpose-plan.md`）。影响：`dot`/`normalize`/`resetall` 里凡触发 dense 共轭的路径在大 n 下仍慢，且受 `kMaxDenseQubits=12` 限制（枚举路径）。

5. **矩阵代数紧凑化**（`MatrixMultiply` 等矩阵×矩阵路径）——LimTDD 侧自认的下一块硬骨头，`Conjugate`/`Transpose` 之外。**已由交叉校验确认是实际 bug**：见 `limtdd-matrixmultiply-bug-report.md`（`MatrixMultiply` 对同尺寸矩阵抛 `size mismatch`（CSX），并对 `ctrl>tgt` 的 S·C·S SWAP 共轭静默算错）。

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

## 5. 2026-08-30 重同步验证记录

**背景**：LimTDD 仓库此前 checkout 在 `temp_debug_dj60`（fidelity/张量网络线），后端头文件 `dd/backend/{DDVector,DDMatrix,DDTypes}.hpp` 只在 `limtdd-backend` 分支上，QReach 的 `QREACH_USE_LIMTDD` 集成因此处于"断线"状态。切回 `limtdd-backend`（HEAD `77936fb`）后重建并验证。

**LimTDD 侧已前进（相对 08-17 记忆/文档）**：level→n 约定已改、门/态/Kron 已紧凑化、`kMaxDenseQubits=12` 仅 guard dense 枚举（VQSS 解锁）、单元测试 44→58。剩余 `Conjugate`/`Transpose` dense + 矩阵代数紧凑化（§3.1.4/5）。

**重建步骤**：
1. LimTDD 切 `limtdd-backend`；`USE_LIMTDD=1 BOOST_PATH=../BOOST/boost_1_81_0 make all`（先清掉旧 `libqreach.so`/`test.o`/`test_qreach` 与 LimTDD 旧 `.o`）。
2. `cp libqreach.so python_pkg/`；`cd python_pkg && USE_LIMTDD=1 BOOST_PATH=../BOOST/boost_1_81_0 ../.venv/bin/python -m invoke build-pybind11`。
3. 编译运行 smoke + workflow 测试。

**验证结果（全绿）**：

| 项 | 结果 |
|---|---|
| `libqreach.so`（LimTDD 核心 3×`.o`） | ✅ 链接成功 |
| `pyqreach` import | ✅ `symbolic_available()=False`（确认 LimTDD） |
| `scripts/limtdd_explicit_smoke.cpp`（语义层 Bell） | ✅ `loc2 dims=[4,1]` |
| `workflow_tests/test_lazy_measurement.py` | ✅ PASSED |
| `workflow_tests/test_ts_structure.py` | ✅ 7/7 PASSED |
| `workflow_tests/test_simulation_grover.py`（Grover_5 vs Qiskit） | ✅ dims=(32,1)，主态 `|11101>` 匹配 |
| `workflow_tests/test_simulation_qft.py`（QFT_5 vs Qiskit） | ✅ 32 态均匀，`satisfy(|+++++>)=True` |
| `examples/test_vqss_correct.py` | ✅ **True**（此前被 kMaxDenseQubits 阻塞） |
| `examples/test_vqss_buggy.py` | ✅ 找到反例（`*** VIOLATION ***`） |
| `examples/test_RUS.py`（while_loop+reset+measure） | ✅ 跑通（InnerProduct 稀疏点积修复有效） |
| `examples/test_rus_wp.py`（backward reachability） | ✅ 跑通（GetLevel 约定修复有效） |

**发现的三处 stale 脚本/注释（非 LimTDD 语义问题，待清）**：
- `scripts/dd_backend_smoke.cpp`：仍用旧 level 约定（`MkBasisVector(1,3)` 现在按 n=1 解释 → `MkBasisVector: index out of range` 崩溃）。需改为 n-based 或标记废弃。
- `examples/test_grover.py`：`from qasm_workflow_runner import ...` 导入失败——reorg 后该模块在 `workflow_tests/`，是 stale import，与后端无关。
- `scripts/limtdd_explicit_smoke.cpp` 内 NOTE 注释已过时（"InnerProduct 对正交态抛异常"已修复）。

**注意**：本次构建使 `python_pkg/libqreach.so` 与根 `libqreach.so` 都变成 LimTDD 版本；LimTDD 的 `.o` 仍写入外部仓库 `../LimTDDexpr/LimTDD/DDPackage/dd/`（既有污染）。切回 CFLOBDD 需 `rm libqreach.so && make all` + 重建 pybind11。

**Statevector 稠密交叉校验 oracle（契约 Phase 6 item 2）已建成**：`python_pkg/workflow_tests/test_backend_crosscheck.py`，23 个小酉电路（覆盖 H/X/Y/Z/S/T、CX/CZ/CP/CSX/SWAP/iSWAP/CCX、U3/ry/rz 浮点参数、非相邻门、非平凡初态）× Qiskit `Statevector` 真值，后端无关。LimTDD 下 21 过 + 2 XFAIL（`csx`、`nonadj_cx`，均为 `MatrixMultiply` bug）；CFLOBDD 下应全绿（待重建验证）。

---

## 6. 相关文档

- `backend-replacement-api-contract.md` — Document A（API 契约，LimTDD 侧实现依据）。
- `limtdd-backend-qreach-plan.md` — QReach 侧计划（Phase 0–6 定义，本文档记录其执行结果）。
- `limtdd-backend-implementation-status.md` — LimTDD 侧实现现状（08-13）。
- `limtdd-backend-implementation-brief.md` — 交给 LimTDD agent 的一页简报。
- `limtdd-innerproduct-bug-report.md` — InnerProduct scaling bug 详情与 resolution。
- `limtdd-matrixmultiply-bug-report.md` — `MatrixMultiply`(matrix×matrix) bug 详情（CSX + ctrl>tgt 两 qubit 门），由交叉校验发现。
- `cflobdd-transpose-dag-corruption.md`、`cflobdd-level8-transpose-bug.md` — CFLOBDD transpose/level≥7 bug 背景（dot/normalize fallback 的由来）。
