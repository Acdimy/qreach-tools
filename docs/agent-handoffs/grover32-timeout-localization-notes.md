# Grover32 Timeout Localization Notes

## Commands Run

- `cd /Users/ftdac/thu/qreach-tools && QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ./.venv/bin/python - <<'PY' ... Grover64 parse_qiskit_cir_lazy ... PY`
- `cd /Users/ftdac/thu/qreach-tools && ./.venv/bin/python - <<'PY' ... subprocess timeout=20 for single-it-grover32-plus.qasm ... PY`
- `cd /Users/ftdac/thu/qreach-tools/python_pkg && QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py --qasm benchmark/converted_qasm/single-it-grover32-plus.qasm --max-instructions 180 --timeout-seconds 20 --skip-fixed-post --bisect`
- `cd /Users/ftdac/thu/qreach-tools/python_pkg && QREACH_PARSE_PROFILE=1 QREACH_PARSE_PROFILE_THRESHOLD=0.05 ../.venv/bin/python workflow_tests/debug_grover32_ccx_pathology.py --qasm benchmark/converted_qasm/single-it-grover32-zero.qasm --max-instructions 180 --timeout-seconds 20 --skip-fixed-post --bisect`
- `cd /Users/ftdac/thu/qreach-tools && ./.venv/bin/python - <<'PY' ... print instructions 96..104 for plus/zero QASM ... PY`

## Grover64 Reference

- Result: completed quickly with no profile lines above threshold. Output: `parsed 959 1 0.2544748749351129`.

## Grover32 Plus

- Short timeout result: emitted profile lines then `TIMEOUT` at 20 seconds.
- Last completed profile line: `[parse-profile] idx=100 op=x qubits=[27] locations_in=1 locations_out=1 elapsed=3.247749s total_locations=102`.
- Last completed profile line in bisect timeout cases was likewise idx=100 `op=x qubits=[27]`, with elapsed around 3.24-3.27s and total_locations=102.
- Next QASM instruction after the last completed line: idx=101 `x [26]`; adjacent following instructions are idx=102 `x [28]`, idx=103 `x [29]`, idx=104 `ccx [28, 29, 46]`.
- Smallest failing prefix: `102`.
- Largest adjacent passing prefix: `status=ok prefix_len=101 num_qubits=63 num_gates=101 num_locations=102 num_result_locations=1 time_load_qasm=0.0029213749803602695 time_build_prefix=0.000680458964779973 time_parse=19.345969249960035 time_fixed_post= time_total=19.364509583916515`.
- Smallest failing result: `status=timeout prefix_len=102 num_qubits= num_gates= num_locations= num_result_locations= time_load_qasm= time_build_prefix= time_parse= time_fixed_post= time_total=20.0`.
- Evidence pattern: plus parse time grows sharply through repeated operations on high-index qubits: idx=89 `ccx [22,23,43]` ~0.66s, idx=94 `ccx [24,25,44]` ~1.45s, idx=99 `ccx [26,27,45]` ~3.35s, and intervening single-qubit `x` gates also become multi-second. Prefix 101 completes just under 20s, while adding prefix instruction 101 (`x [26]`) crosses the 20s timeout.

## Grover32 Zero

- Last completed profile line: none; no `[parse-profile]` line exceeded threshold in the bisect run.
- Smallest failing prefix: `None` within `--max-instructions 180`.
- Largest adjacent passing prefix: `status=ok prefix_len=179 num_qubits=63 num_gates=179 num_locations=180 num_result_locations=1 time_load_qasm=0.0026431670412421227 time_build_prefix=0.0012200410710647702 time_parse=0.08474900003056973 time_fixed_post= time_total=0.10402729199267924`.
- Nearby instruction check for zero showed indices 96..104 are CCX-heavy (`ccx [60,61,62]`, `ccx [58,59,61]`, ..., `ccx [44,45,54]`) but still parse rapidly under lazy parsing.

## Current Hypothesis

Grover32 plus timeout is not caused by CCX count alone and does not reproduce on the zero variant. The slowdown begins after the plus variant's fan-in/tree setup creates complex lazy quantum operations on higher qubit pairs; subsequent `x` operations on those already-complex qubits become as expensive as the preceding `ccx`. The smallest 20s timeout prefix is 102, but the first clearly pathological growth is already visible around idx=89/94/99, culminating after idx=100 with the next `x [26]` crossing the timeout.

## Next Inspection Target

Lazy post-image / `QOperation` gate application path for single-qubit X and CCX on accumulated plus-state Grover32 terms, then C++ CFLOBDD gate construction/application if Python only dispatches directly.
