We want to extend QisMC's existing counterexample output with a lightweight
semantic debugging diagnosis mechanism. Please first inspect the existing
implementation carefully, identify the relevant SP/WP fixed-point code,
counterexample-generation code, transition-system representation, and
location/source mapping, and then implement the feature with minimal changes
to the current verification architecture.

Do NOT redesign the qCTL semantics, SP/WP definitions, or the existing
bidirectional fixed-point algorithm unless there is a genuine implementation
bug. The goal is to reuse information already available during verification.

======================================================================
1. Current behavior
======================================================================

For the buggy RUS motivating example, QisMC currently produces a
counterexample similar to:

=== Counterexample Analysis ===
Violation: location 18 (identifier: 'S4.EW') satisfies 'outloop' but NOT 's'.

Trace (location → identifier → edge):
 Step  1: location 0   id='S0'       --X(2)[]-->                 [Zero=False outloop=False s=False]
 Step  2: location 1   id='S1'       --meas1(2)[]-->
 Step  3: location 3   id='S2'       --meas0(1)[]reset0(1)[]-->
 Step  4: location 5   id='S3'       --I-->
 Step  5: location 7   id='S4.W'     --H(0)[]-->
 Step  6: location 8   id='S4.W.S1'  --T(0)[]-->
 Step  7: location 9   id='S4.W.S2'  --CX(0,1)[]-->
 Step  8: location 10  id='S4.W.S3'  --H(0)[]-->
 Step  9: location 11  id='S4.W.S4'  --CX(0,1)[]-->
 Step 10: location 12  id='S4.W.S5'  --T(0)[]-->
 Step 11: location 13  id='S4.W.S6'  --H(0)[]-->
 Step 12: location 14  id='S4.W.S7'  --meas0(0)[]-->
 Step 13: location 15  id='S4.W.S8'  --I-->
 Step 14: location 17  id='S4.WN'    --I-->
 Step 15: location 18  id='S4.EW'    ---->                      [outloop=True] *** VIOLATION ***

This output is useful as a violating execution witness, but it does NOT
automatically identify the root cause. In particular, the specification does
not explicitly say "there should be a reset at the loop entry", so QisMC
must NOT claim that the counterexample itself proves that the missing reset
is the bug.

We want to provide additional semantic information that helps the programmer
localize a suspicious program region.

======================================================================
2. RUS setup for the new diagnosis
======================================================================

For this motivating example, use TWO quantum annotations:

Initial annotation:
    f(l_init) = Span(|000>)

Exit/target annotation:
    f(l_exit) = Span(
        1/sqrt(3) |001> + i*sqrt(2)/sqrt(3) |011>
    )

Use the existing bidirectional fixed-point machinery to compute the ordinary:

    sp[l]   strongest postcondition / reachable subspace
    wp[l]   weakest precondition / backward admissible subspace

Do NOT introduce a new "diagnostic WP", "DWP", or a new semantic object.
We intentionally want to keep the theory and implementation simple.

The final fixed-point WP values should be exactly the WP values already
defined by QisMC's existing annotation semantics.

======================================================================
3. Important observation: final SP/WP mismatch is NOT sufficient
======================================================================

Do NOT simply search the final fixed point for:

    sp[l] ⊄ wp[l]

and report the first such location as the bug location.

This is incorrect as a localization strategy.

Once an incompatibility exists downstream, WP propagation can cause many
earlier locations (potentially a large portion of the loop) to satisfy

    sp[l] ⊄ wp[l].

Therefore:

    "first inconsistent location in the final SP/WP result"

is NOT equivalent to:

    "location where the bug is introduced".

The new diagnosis must instead inspect HOW SP is constructed during the
forward fixed-point computation.

======================================================================
4. Core idea: SP update provenance
======================================================================

After verification returns False:

1. Keep/freeze the FINAL wp[l] computed by the ordinary bidirectional
   fixed-point algorithm.

2. Inspect or replay the forward construction of sp[l].

For every forward propagation/update caused by an incoming transition

    u = (p --E--> l)

let:

    C = SP_E(sp[p])

be the subspace contribution propagated through this edge.

Before this update, let:

    S_old = sp[l]

After joining the contribution:

    S_new = S_old ∨ C

where ∨ is the existing subspace join operation.

We are interested in updates satisfying:

    S_old ⊆ wp[l]

but

    S_new ⊄ wp[l].

In words:

    Before this propagation event, all states accumulated at location l
    were compatible with the backward-computed requirement.

    This particular update expands the reachable subspace and introduces
    states that are no longer fully contained in wp[l].

Call this internally something simple such as:

    compatibility-breaking SP update

or

    incompatible SP expansion

Do not overstate it as a mathematically proven "root cause".

======================================================================
5. Expected behavior on the RUS example
======================================================================

The buggy RUS program misses a reset of the measured ancilla qubit when the
execution starts another loop iteration.

Semantically, the important behavior is:

- On the first entry into the loop, the reachable subspace at the loop-entry
  location is one-dimensional.

- After a failed RUS iteration, the measurement leaves a residual state that
  travels along the loop back-edge.

- Because the reset is missing, this residual state is not reinitialized
  before the next iteration.

- When the loop-entry location is updated for the second time through the
  back-edge, its SP expands, expectedly from dimension 1 to dimension 2.

The diagnosis we expect is approximately:

    Location: S4.W (loop entry)
    Incoming contribution: loop back-edge
    SP dimension: 1 -> 2
    Before update:
        S_old ⊆ wp[S4.W]       true
    After update:
        S_new ⊆ wp[S4.W]       false

This should identify the loop re-entry / back-edge propagation as a
semantically suspicious event.

This DOES NOT mean QisMC has automatically inferred:

    "insert reset(q0) here"

Instead, QisMC should tell the programmer:

    this is the first observed SP expansion that introduces reachable
    quantum behavior incompatible with the future annotation.

The programmer can then inspect the corresponding source region and discover
the missing reset.

======================================================================
5.5. Semantic alignment of the current RUS TS with the paper-level story
======================================================================

While inspecting the current implementation, be careful about what the
parser-generated location names mean.

In the present explicit transition-system encoding produced by
`python_pkg/parse_qiskit.py`, the location named:

        S4.W

is the control-flow entry into the while-body *before* the first instruction
inside that body is executed.

For the corrected RUS circuit, the first instruction inside the loop body is:

        reset(0)

Therefore, in the current implementation:

- `S4.W` is a pre-reset loop-entry location;
- the first post-reset location (currently `S4.W.S1` in the corrected
    encoding) is the better semantic match for the paper-level intuition
    "the loop starts a new iteration with the ancilla reinitialized".

This matters for interpreting dimensions:

- in the corrected circuit, the current TS may still allow the reachable
    subspace at `S4.W` itself to contain residual states arriving from the
    previous failed iteration;
- the reset can then collapse that subspace back to the intended
    one-dimensional reinitialized state at the immediate successor location;
- in the buggy circuit, the corresponding reset edge is absent, so this
    reinitialization point disappears and the incompatible residual component is
    allowed to continue propagating through the loop body.

Consequently, if the debugging output is compared against the paper-level
statement

        "first loop entry is 1-dimensional; after the back-edge it becomes 2-dimensional"

the implementation should NOT assume that the raw parser identifier `S4.W`
must be that exact semantic point.

Instead, the diagnosis should align to the nearest semantic point consistent
with the implementation:

- if the loop body begins with a reset/initialization, the semantically
    reinitialized iteration point is the first post-reset location, not the
    pre-reset control-flow entry;
- if the loop body does not begin with a reset/initialization, the raw loop
    entry itself remains the relevant comparison point.

For the initial implementation, it is acceptable to report the actual
implementation locations explicitly, for example:

- pre-reset loop entry: `S4.W`
- first post-reset location in the corrected program: `S4.W.S1`

and to explain in the diagnostic text which one is being compared to the
paper-level intuition.

Do not treat this alignment issue as a semantic contradiction in the SP/WP
theory.  It is primarily a modelling / source-mapping issue in the current TS
naming scheme, and it should be documented as such.

======================================================================
6. Worklist-order issue
======================================================================

Be careful about the phrase "first".

The SP fixed-point algorithm may use a worklist, and the exact order in which
updates are processed may affect which update is encountered first.

Please inspect the current implementation and determine:

- Is the worklist/update order deterministic?
- Does the RUS example consistently expose the loop-back-edge update described
  above?
- Can we record the provenance of each SP expansion reliably?
- Can we associate each expansion with the incoming transition that caused it?

For the initial implementation, deterministic algorithmic order is acceptable,
but document this clearly in the code.

If straightforward, also compute an edge-level diagnostic:

    C_e = SP_E(sp[p])

and test:

    C_e ⊆ wp[l]

for incoming contributions.

This can help distinguish which predecessor/back-edge contributes incompatible
states at a join point and may be more stable than relying solely on worklist
ordering.

However, DO NOT substantially redesign the fixed-point algorithm just to make
this canonical. Keep the implementation lightweight.

======================================================================
7. Information to record
======================================================================

Please determine the cheapest way to record enough SP provenance during
verification, or replay SP propagation after verification.

For each relevant SP-changing update, ideally retain:

- target location index;
- target location identifier;
- predecessor location;
- incoming edge / quantum operation;
- iteration/update sequence number;
- dimension of S_old;
- dimension of propagated contribution C;
- dimension of S_new;
- whether S_old ⊆ wp[l];
- whether S_new ⊆ wp[l];
- if feasible, whether C ⊆ wp[l];
- whether the update is caused by a loop back-edge / CFG back-edge;
- source-level mapping if already available.

Do NOT print full exponentially large matrices/state vectors by default.

Dimension, containment relationships, edge identity, and source locations are
the primary debugging information.

If symbolic subspaces have a concise representation in the current backend,
consider making detailed subspace output optional under a verbose/debug flag.

======================================================================
8. Desired output
======================================================================

Keep the existing counterexample output.

Add a separate section after it, for example:

=== Semantic Diagnosis ===

Potential incompatibility-introducing propagation:

  Location:
      7 ('S4.W')

  Incoming edge:
      <loop back-edge / predecessor information>

  SP update:
      dimension 1 -> 2

  Compatibility:
      before update: SP ⊆ WP   [True]
      after update:  SP ⊆ WP   [False]

  Incoming contribution:
      contribution ⊆ WP        [False]

  Interpretation:
      This propagation introduces reachable quantum states that are
      incompatible with the backward requirement at the loop entry.
      Inspect the corresponding loop re-entry path and source region.

Avoid output such as:

    "Bug found at S4.W"
    "Missing reset detected"
    "Root cause: reset"

unless such information is independently justified by another analysis,
which currently it is not.

======================================================================
9. Architecture constraints
======================================================================

The implementation should be minimally invasive.

Prefer something like:

existing verification
    |
    +-- final SP/WP fixed point
    |
    +-- CTL/qCTL result
    |
    +-- NuSMV counterexample
    |
    +-- if result == False:
            semantic diagnosis
                |
                +-- inspect/replay SP update history
                +-- compare updates with final WP
                +-- rank/report suspicious propagation events

Do NOT modify:

- qCTL syntax;
- qCTL semantics;
- definition of quantum atomic propositions;
- mathematical definitions of SP/WP;
- qCTL-to-CTL reduction;
- NuSMV semantics;
- existing correctness assumptions.

Reuse the existing implementations of:

- SP transformers;
- WP transformers;
- measurement semantics;
- initialization/reset semantics;
- unitary transformations;
- subspace join/meet;
- subspace containment/equality;
- dimension computation.

======================================================================
10. Scope and limitations
======================================================================

This diagnosis is intentionally specification-dependent.

It works particularly well when the program has both:

- a meaningful initial quantum annotation, which constrains forward SP; and
- a meaningful later/exit quantum annotation, which constrains backward WP.

If the future annotation is absent and WP remains the entire Hilbert space
(Identity), then:

    sp[l] ⊆ wp[l]

may always hold and this diagnosis may provide little or no localization
information.

This is acceptable.

Do NOT introduce new semantics just to handle this case.

The current research goal is NOT a universal quantum fault-localization
algorithm. It is a lightweight counterexample-guided semantic diagnosis
mechanism that provides richer debugging information when the specification
contains sufficient quantum constraints.

======================================================================
11. What I want you to do first
======================================================================

Before changing code:

1. Locate the implementation of the bidirectional SP/WP fixed-point algorithm.
2. Explain how SP updates are currently scheduled and joined.
3. Locate the implementation of:
   - measurement SP/WP;
   - initialize/reset SP/WP;
   - unitary SP/WP;
   - subspace containment;
   - subspace dimension.
4. Locate the counterexample analysis/output code.
5. Determine whether SP update provenance can be recorded cheaply, or whether
   replaying the SP fixed-point after verification is cleaner.
6. Trace the buggy RUS example through the CURRENT implementation and verify
   empirically whether:
       - first loop entry has SP dimension 1;
       - the loop back-edge later expands it to dimension 2;
       - the pre-update SP is contained in final WP at that location;
       - the post-update SP is not contained in final WP.
7. Verify the above with actual symbolic subspaces, not only dimensions.
   Dimension growth alone is NOT sufficient evidence of incompatibility.
8. Report your findings before making large architectural changes.

Then propose the smallest implementation plan.

After that, implement the diagnosis, add tests using the RUS example, and show
the resulting output.

======================================================================
12. Testing requirements
======================================================================

At minimum add/check tests for:

A. Buggy RUS with both annotations:
   Expected to produce a compatibility-breaking update near loop re-entry.

B. Correct RUS with the reset restored:
   The same suspicious propagation should disappear, or the relevant SP
   contribution should remain compatible with WP.

C. RUS without the final target annotation:
   Diagnosis should gracefully report that no useful SP/WP incompatibility was
   found, rather than failing.

D. Existing benchmarks:
   Verification results must remain unchanged by enabling/disabling diagnosis.

The diagnosis layer must not change model-checking truth values.

======================================================================
13. Research intent
======================================================================

This implementation is intended to support a paper-level claim approximately
like:

"After detecting a property violation, QisMC further exploits its
bidirectional analysis for debugging. It fixes the backward-computed weakest
preconditions and examines the incremental construction of strongest
postconditions. When an SP update expands the reachable subspace beyond the
corresponding WP, QisMC highlights the responsible propagation event and its
source-level control-flow region. This does not in general identify the root
cause automatically, but provides semantic guidance for inspecting the
violating execution."

Please keep this intended claim in mind when designing the implementation:
the code should provide enough concrete evidence to support this claim, but
should not implement unnecessary generalizations that make the existing QisMC
architecture or theory substantially more complicated.