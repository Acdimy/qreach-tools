# QisMC Weakest-Precondition Computation

This note summarizes the weakest-precondition (WP) computation defined in
the QisMC paper. The new debugging diagnosis should reuse this existing
semantics directly. Do NOT introduce DWP or another backward semantics.

## 1. Meaning of WP

For each program location \(l\), QisMC computes a subspace \(wp[l]\).

Intuitively:

- \(sp[l]\) describes the quantum states that can reach \(l\) from preceding
  annotations.
- \(wp[l]\) describes the quantum states at \(l\) that are compatible with
  subsequent annotations.

Thus SP propagates information forward, while WP propagates annotation
constraints backward.

## 2. Backward fixed-point computation

Let \(X=(X_l)_{l\in L}\) be the current collection of subspaces.

QisMC defines the backward update at location \(l\) as

\[
Y_l
=
X_l
\wedge
\left(
\bigwedge_{k\in succ(l)}
\mathcal E_{lk}^{-1}(X_k)
\right)
\wedge
\bigwedge \widetilde f(l),
\]

where:

- \(succ(l)\) is the set of successors of \(l\);
- \(\mathcal E_{lk}\) is the quantum operation on transition \(l\rightarrow k\);
- \(\mathcal E_{lk}^{-1}(X_k)\) is the backward pre-image of \(X_k\);
- \(\widetilde f(l)\) is the set of quantum annotations at \(l\);
- \(\wedge\) denotes subspace intersection.

If no quantum annotation is attached to \(l\),

\[
\bigwedge\widetilde f(l)=\operatorname{supp}(I)=\mathcal H.
\]

WP is obtained through a greatest fixed-point iteration. Consequently, the
subspaces monotonically shrink:

\[
wp^{(0)}[l]
\supseteq
wp^{(1)}[l]
\supseteq
\cdots
\supseteq
wp^*[l].
\]

## 3. Unitary transition

If \(\mathcal E_{lk}\) is a unitary operation \(U\), its pre-image is the
ordinary inverse transformation:

\[
\mathcal E_{lk}^{-1}(X)=U^\dagger X.
\]

Equivalently, for

\[
X=\operatorname{span}\{
|\phi_1\rangle,\ldots,|\phi_m\rangle
\},
\]

we have

\[
U^{-1}(X)
=
\operatorname{span}\{
U^\dagger|\phi_1\rangle,\ldots,
U^\dagger|\phi_m\rangle
\}.
\]

A unitary therefore preserves the dimension of the subspace during backward
propagation.

## 4. Projective measurement branch

If \(\mathcal E_{lk}\) corresponds to a measurement branch represented by
projector \(P\), QisMC defines

\[
\boxed{
P^{-1}(X)
=
(P\wedge X)\vee(I-P)
}
\]

where \(\vee\) denotes subspace join/span.

The intuition is:

- the component inside \(P\) takes this measurement branch and therefore
  must satisfy \(X\);
- the component inside \(I-P\) does not take this branch and consequently
  imposes no constraint for this particular transition.

For a computational-basis measurement with projectors \(P_0\) and \(P_1\),

\[
P_0^{-1}(X)
=
(P_0\wedge X)\vee P_1,
\]

and

\[
P_1^{-1}(X)
=
(P_1\wedge X)\vee P_0.
\]

If a measurement produces multiple outgoing transitions, the backward
constraints from all successors are combined by intersection in the global
WP update.

## 5. Initialization

If \(\mathcal E_{lk}\) initializes the relevant quantum register to a fixed
state \(|\psi\rangle\), QisMC defines

\[
\boxed{
\mathcal E_{lk}^{-1}(X)
=
\begin{cases}
0, & |\psi\rangle\notin X,\\
I, & |\psi\rangle\in X.
\end{cases}
}
\]

The intuition is that initialization discards the previous state and prepares
\(|\psi\rangle\).

Therefore:

- if the initialized state satisfies the required postcondition, every
  possible input is admissible;
- otherwise, no input can satisfy the required postcondition.

A reset is the special case of initialization to \(|0\rangle\).

When initialization/reset acts only on a subsystem, follow the existing QisMC
implementation for embedding this semantics into the complete Hilbert space.
Do not reinterpret the paper's notation as replacing the entire global state.

## 6. Multiple successors

For

\[
l\rightarrow k_1,\ldots,l\rightarrow k_m,
\]

each successor contributes a backward constraint

\[
C_i
=
\mathcal E_{lk_i}^{-1}(wp[k_i]).
\]

These constraints are combined by intersection:

\[
wp[l]
\leftarrow
wp[l]
\wedge C_1
\wedge\cdots\wedge C_m.
\]

This gives the important SP/WP duality:

\[
\boxed{
SP:\quad\text{forward image + join over reachable behaviors}
}
\]

\[
\boxed{
WP:\quad\text{backward pre-image + intersection over future constraints}
}
\]

## 7. Use in the new debugging diagnosis

The debugging extension must NOT modify the WP computation above.

First compute the ordinary final fixed-point result

\[
wp^*[l].
\]

Then inspect or replay the incremental construction of SP.

Suppose an incoming transition

\[
p\xrightarrow{\mathcal E}l
\]

produces the forward contribution

\[
C=\mathcal E(sp[p]).
\]

For an SP update, define

\[
S_{\mathrm{old}}=sp[l],
\]

and

\[
S_{\mathrm{new}}
=
S_{\mathrm{old}}\vee C.
\]

The diagnostic event of interest satisfies

\[
\boxed{
S_{\mathrm{old}}\subseteq wp^*[l]
\quad\land\quad
S_{\mathrm{new}}\nsubseteq wp^*[l].
}
\]

This means that the incoming propagation introduces reachable quantum states
that are incompatible with the future constraints captured by the existing
WP.

This is a debugging signal, not a new WP semantics and not necessarily the
root cause of the bug.

## 8. RUS interpretation

For the RUS motivating example, we use two annotations:

\[
\widetilde f(l_{\mathrm{init}})
=
\operatorname{Span}(|000\rangle),
\]

and

\[
\widetilde f(l_{\mathrm{exit}})
=
\operatorname{Span}
\left(
\frac{1}{\sqrt3}|001\rangle
+
\frac{i\sqrt2}{\sqrt3}|011\rangle
\right).
\]

The initial annotation provides meaningful forward SP information, while the
target annotation provides meaningful backward WP information.

The expected buggy behavior is that the first loop-entry SP remains compatible
with its WP, whereas the contribution propagated through the loop back-edge
after a failed iteration expands the reachable subspace beyond that WP.

Implementation alignment note for the current explicit TS:

- the parser-generated location `S4.W` is the control-flow entry into the
  while-body before the first statement in that body executes;
- in the corrected RUS program, the first statement is `reset(0)`, so the
  first post-reset location is the better semantic match for the paper-level
  phrase "re-enter the loop with the ancilla reinitialized";
- therefore, a current implementation trace may show `S4.W` itself carrying a
  residual higher-dimensional subspace in the corrected program, while the
  immediate post-reset successor collapses back to the intended one-dimensional
  state.

This should be treated as a modelling/source-mapping alignment issue, not as
an intended difference in SP/WP semantics.  When comparing implementation
results against the paper-level RUS explanation, align the diagnosis to the
nearest semantically equivalent program point rather than assuming the raw
location name `S4.W` is already that point.

The diagnosis should therefore inspect actual subspace containment:

\[
S_{\mathrm{old}}\subseteq wp[l],
\qquad
S_{\mathrm{new}}\nsubseteq wp[l].
\]

A dimension change such as

\[
\dim S_{\mathrm{old}}=1
\rightarrow
\dim S_{\mathrm{new}}=2
\]

is useful diagnostic information, but dimension growth alone is NOT evidence
of incompatibility.

## 9. Implementation principle

The diagnosis implementation should:

1. reuse the existing final WP fixed point;
2. reuse the existing unitary, projector, and initialization pre-images;
3. preserve the existing SP/WP fixed-point algorithms;
4. record or replay SP propagation provenance;
5. compare each relevant SP expansion against the corresponding final WP;
6. report the incoming edge/location responsible for a compatibility-breaking
   expansion;
7. leave all original model-checking results unchanged.

If the current implementation appears inconsistent with the formulas above,
report the discrepancy before modifying the verification semantics.