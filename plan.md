# qreach-tools Current Plan

## Status

The original near-term goal of getting a working SymTS post-image pipeline is complete.

Completed items include:

- classical proposition support inside `transition_system_qadd.hpp`
- Python exposure through `qreach_python_wrapper.cpp`
- minimal symbolic-vs-explicit post regression
- reduced control-flow validation for `while_loop` and `if_test`
- Steane-derived bounded regression harness
- fused symbolic post prototype that removes explicit `exists_vars` and `rename_vars` from the hot path

The remaining work is no longer about making the symbolic path merely functional. It is now about scaling, completing the higher-level model-checking workflow, and turning current profiling results into structural optimizations.

## Remaining Work From The Original Roadmap

### 1. Complete the full SymTS qCTL/model-checking workflow

Goal:

- connect the current symbolic transition system to the existing qCTL labeling and checking workflow in `python_pkg/qctl.py`
- reproduce something closer to the original `test_steane_qiskit.py` end-to-end behavior using SymTS, not just symbolic post fixed points

Why it still matters:

- the repository already has the logical labeling layer, but the current symbolic work only covers reachability/post propagation
- this is the key functional milestone still missing from the original broader plan

### 2. Add stronger symbolic-vs-explicit regression coverage

Goal:

- extend current regressions beyond minimal and bounded Steane checks
- add small but semantically rich symbolic-vs-explicit checks for additional measurement/control-flow patterns

Why it still matters:

- current bounded Steane regression protects structure and stability, but it does not yet cover full qCTL semantics or larger symbolic equivalence surfaces

### 3. Decide whether to keep the current fused post as the long-term kernel

Goal:

- treat the current fused post implementation as the baseline kernel for further optimization
- avoid going back to the old `Apply -> exists -> rename` path unless a correctness issue appears

Why it still matters:

- the fused post path clearly wins on small and medium workloads by removing the old abstraction and renaming hotspot
- the next design work should build on this kernel, not branch the implementation surface again

## Optimization Directions That Still Make Sense

### A. Build a source-indexed relation-side post kernel

Goal:

- replace or augment the current fused recursion with a structure indexed by source-location blocks or source prefixes
- avoid re-traversing the full interleaved relation DD for every delta fragment and for every fixed-point round

Expected benefit:

- attack the real remaining bottleneck directly, since the present Steane `measure_count=1` timings are dominated by recursive `Apply(APPLY)` structure traversal

Why this is the highest-value next optimization:

- recent profiling already ruled out `exists`, `rename`, terminal `postImage`, and several local memo/short-circuit ideas as the dominant remaining issue

### B. Add a precompiled relation cache keyed by source encoding prefixes

Goal:

- precompute a relation-side structure that groups outgoing transitions by source-bit prefixes or by exact source locations
- reduce repeated work across the 38 to 40 fixed-point iterations seen in Steane `measure_count=1`

Expected benefit:

- turn the current repeated recursive scanning of the interleaved relation into more direct outgoing-transition retrieval

Risk:

- this is more invasive than local `Apply` short-circuits and may increase memory usage

### C. Add profiling counters that expose relation recursion width, not only time

Goal:

- record counts such as fused recursive calls, memo hits, zero-subtree prunes, and target-node reconstructions

Expected benefit:

- make it easier to distinguish “too many recursive calls” from “too expensive recursive calls”
- prevent future optimization work from relying only on wall-clock timing

### D. Add regression coverage for fused post structural invariants

Goal:

- keep the current bounded Steane regression
- add regression assertions for iteration counts or bounded-node growth where stable enough

Expected benefit:

- protect future refactors of the symbolic kernel from silently regressing the current fused-post behavior

## Future Functional Work

### Full symbolic pre-image or bidirectional model checking

Goal:

- decide whether symbolic pre-image is worth implementing after the post-only path stabilizes

Reason:

- the original explicit workflow has both pre and post model-checking capabilities, while the symbolic path is still post-centric

### Better SymTS integration into existing Python workflows

Goal:

- make `parse_qiskit.py`, `qctl.py`, and the symbolic wrapper feel like one workflow instead of several stitched pieces

Reason:

- this reduces friction for future benchmarks and avoids special-case test harnesses for every symbolic experiment

### Cleanup and hardening

Goal:

- remove dead experimental branches once the preferred symbolic kernel is clear
- document profiling knobs and regression entrypoints in a more durable way

Reason:

- the repository now contains several generations of symbolic-post experiments and bounded debugging scripts

## Immediate Recommendation

If development continues from the current version, the best next engineering step is:

- design a source-indexed or source-grouped specialized relation-side post kernel for `transition_system_qadd.hpp`

This is the most consistent conclusion from the current measurements. The remaining runtime is dominated by recursive traversal of the interleaved relation structure across many fixed-point iterations, not by terminal operations or the previously suspected abstraction/renaming phases.