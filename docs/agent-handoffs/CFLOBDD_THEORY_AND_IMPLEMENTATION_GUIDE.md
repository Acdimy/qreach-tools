# CFLOBDD Theory and Implementation Guide for AI Coding Agents

> Purpose: This document is a repository-oriented guide for reading, modifying, debugging, and extending a CFLOBDD implementation without requiring direct access to the original paper.
>
> Scope: theory, representation invariants, runtime architecture, canonicalization, interpretation, construction, Apply/PairProduct/Reduce, matrix operations, and implementation hazards.
>
> Primary source: *CFLOBDDs: Context-Free-Language Ordered Binary Decision Diagrams* (Sistla, Chaudhuri, and Reps).

---

## 1. Executive Summary

A CFLOBDD is a canonical hierarchical decision diagram for functions

\[
f:\{0,1\}^{2^k}\to V
\]

where `k` is the CFLOBDD level and `V` is a finite terminal-value domain.

It is intended to be a plug-compatible alternative to BDDs/MTBDDs/ADDs, but its internal sharing principle is different:

- A BDD shares suffix sub-DAGs.
- A CFLOBDD shares hierarchical components that behave like non-recursive procedures.
- The same lower-level component can be called from multiple contexts.
- Call and return edges must match, so legal paths form a balanced-parenthesis context-free language.
- A level-`k` grouping consumes exactly `2^k` Boolean variables along every legal path.
- Variables are interpreted contextually: the `i`-th visit to a level-0 grouping reads the `i`-th variable in the fixed ordering.
- Canonical structural invariants plus hash-consing allow semantic equality to become pointer equality.

The central implementation pattern is:

```text
CFLOBDD
  = root Grouping
  + ValueTuple

Internal Grouping
  = AConnection
  + AReturnTuple
  + array of BConnections
  + array of BReturnTuples
  + numberOfExits
```

Operationally:

```text
interpret(grouping, assignment):
    split assignment into first half and second half
    i = interpret(AConnection, first half)
    j = interpret(BConnections[i], second half)
    return BReturnTuples[i][j]
```

The most important nontrivial operation is binary Apply:

```text
PairProduct
    -> compute terminal result tuple
    -> collapse duplicate values left-to-right
    -> Reduce the product grouping according to that collapse
    -> hash-cons the final result
```

Do not treat CFLOBDD reduction as BDD-style local node elimination. CFLOBDD reduction is a recursive backwards propagation of exit equivalence classes through B-connections and then the A-connection.

---

## 2. Mental Model

### 2.1 BDD Mental Model

A BDD is an acyclic branching program:

```text
test x0
  -> test x1
      -> ...
          -> terminal
```

Each node is tied to a fixed variable level. Sharing occurs when different prefixes reach the same suffix computation.

### 2.2 CFLOBDD Mental Model

A CFLOBDD is better understood as a single-entry, multi-exit, non-recursive hierarchical finite-state machine.

A level-`k` grouping behaves like a procedure:

```text
procedure G_k(input bits[0 .. 2^k-1]) -> exit index:
    i = call A on first 2^(k-1) bits
    j = call B_i on last 2^(k-1) bits
    return BReturnTuple_i[j]
```

The A-call decides which middle vertex is reached. Each middle vertex selects one B-call. The B-call decides which exit of the current grouping is reached.

The graph may look cyclic if call/return edge types are ignored. Semantically it is not cyclic because only properly matched call/return paths are legal.

### 2.3 Why the Compression Can Be Stronger Than BDDs

Suppose a recursively defined object repeatedly uses the same lower-level function twice, e.g., a Kronecker-recursive matrix.

A BDD may need multiple copies or multiple distinct residual subgraphs.

A CFLOBDD can add one constant-size outer grouping that calls the same previous-level grouping:

```text
G_{k+1}:
    AConnection -> G_k
    BConnection[1] -> G_k
    BConnection[2] -> G_k
```

For favorable families, size grows as `O(k)` while the represented decision tree has doubly exponential size in `k`.

This is especially effective for:

- tensor/Kronecker-recursive matrices,
- identity and Hadamard-like structures,
- structured relations,
- quantum state vectors and operators,
- large regular Boolean functions.

It is not a universal guarantee. Arbitrary functions can still produce large CFLOBDDs.

---

## 3. Core Data Types

A typical implementation contains equivalents of the following classes.

```cpp
abstract class Grouping {
    int level;
    int numberOfExits;
};

class InternalGrouping : public Grouping {
    GroupingHandle AConnection;
    ReturnTupleHandle AReturnTuple;

    std::vector<GroupingHandle> BConnections;
    std::vector<ReturnTupleHandle> BReturnTuples;
};

class ForkGrouping : public Grouping {
    level = 0;
    numberOfExits = 2;
};

class DontCareGrouping : public Grouping {
    level = 0;
    numberOfExits = 1;
};

template<class Value>
class CFLOBDD {
    GroupingHandle grouping;
    ValueTuple<Value> valueTuple;
};
```

Repository names may differ. Search for concepts rather than exact names:

| Paper term | Likely code names |
|---|---|
| Grouping | `Node`, `Grouping`, `NodeHandle`, `NWAOBDDNode`, `CFLOBDDNode` |
| InternalGrouping | `InternalNode`, `InternalNodeHandle` |
| ForkGrouping | `ForkNode`, `ForkNodeHandle` |
| DontCareGrouping | `DontCareNode`, `NoDistinctionNode` |
| ReturnTuple | `ConnectionReturnMap`, `ReturnMap`, `IntPair`, `Tuple`, `ReductionMap` |
| ValueTuple | `ValuesList`, `WeightedValuesList`, `MapHandle` |
| RepresentativeGrouping | `Canonicalize`, `MkNode`, `Intern`, `GetCanonical`, `HashCons` |
| NoDistinctionProtoCFLOBDD | `NoDistinctionNode`, `MkNoDistinction`, `GetNoDistinctionNode` |
| PairProduct | `PairProduct`, `CrossProduct`, `ApplyProduct` |
| Reduce | `Reduce`, `ReduceNode`, `ReduceExits` |

### 3.1 Level-0 Groupings

There are exactly two semantic base groupings.

#### Fork grouping

```text
input bit 0 -> exit 1
input bit 1 -> exit 2
```

#### Don't-care grouping

```text
input bit 0 -> exit 1
input bit 1 -> exit 1
```

Important: neither object is permanently attached to a variable name.

---

## 4. Matched Paths

### 4.1 Call/Return Matching

Every connection from a level-`k` grouping to a level-`k-1` grouping is a logical call edge. Its corresponding return edges are tied to that particular call site.

A legal path must obey stack discipline:

```text
call A
    call X
    return X
return A
```

but not:

```text
call A
    call X
return A
return X
```

The legal edge-label strings form a balanced-parenthesis language.

This explains the name:

```text
CFL + OBDD = Context-Free-Language Ordered Binary Decision Diagram
```

### 4.2 Implementation Consequence

Return edges conceptually belong to the caller, not the callee.

The callee only exposes numbered exits. The caller stores a tuple saying where each callee exit returns.

```text
callee exit index
    --caller return tuple-->
caller middle/exit index
```

Therefore a lower-level grouping can be reused in different contexts with different return tuples.

This is the principal structural advantage over BDD suffix sharing.

---

## 5. Contextual Variable Interpretation

A level-0 grouping is not labeled with a fixed variable.

For a fixed global variable ordering:

```text
x0, x1, ..., x_{2^k-1}
```

the `i`-th level-0 visit on a legal path reads `x_i`.

Consequences:

1. The same `ForkGrouping` singleton can represent every variable occurrence.
2. Variable position is encoded by hierarchical call context.
3. Reordering variables is a structural operation, not a relabeling of leaf nodes.
4. Code that assumes a grouping has a permanent variable index is conceptually wrong.

---

## 6. Structural Semantics

### 6.1 Number of Variables Per Level

Let `d(k)` be the number of level-0 decisions encountered on a matched path through a level-`k` grouping.

```text
d(0) = 1
d(k) = d(k-1) + d(k-1)
```

Therefore:

```text
d(k) = 2^k
```

A level-`k` CFLOBDD represents a function of exactly `2^k` Boolean arguments.

### 6.2 Exit Meaning

A grouping does not directly map assignments to terminal values.

Instead it maps assignments to exit indices:

```text
Grouping: assignment -> exit index
```

The top-level CFLOBDD maps exit indices to values:

```text
ValueTuple[exit index] -> Value
```

This separation is critical:

- `Grouping` captures partition structure.
- `ValueTuple` labels partition classes.
- The same grouping can be reused with different value tuples.
- Unary operations often need only transform the value tuple.

### 6.3 Denotational View

A grouping with `m` exits denotes a partition:

```text
L[1], L[2], ..., L[m]
```

where `L[i]` is the set of assignments reaching exit `i`.

The partition sets are disjoint and cover the entire assignment space.

For level 0:

```text
ForkGrouping     = [{0}, {1}]
DontCareGrouping = [{0,1}]
```

At an internal grouping, languages from A and B components are concatenated, then routed through B return tuples.

This view is useful for proving correctness of PairProduct:

```text
product exit (i,j)
    denotes
assignments reaching exit i in left operand
INTERSECT
assignments reaching exit j in right operand
```

---

## 7. Return Tuples

### 7.1 Meaning

A return tuple is an indexed map.

For a B-connection:

```text
BReturnTuple[j]
```

maps exit `j` of the called B grouping to an exit of the current grouping.

For the A-connection:

```text
AReturnTuple[j]
```

maps exit `j` of the called A grouping to a middle vertex of the current grouping.

Typical representation:

```cpp
std::vector<unsigned> returnTuple;
```

Paper indexing is 1-based. C++ implementations are often 0-based. Be extremely careful when porting algorithms.

### 7.2 Value Tuple

The top-level value tuple maps root exits to values:

```cpp
valueTuple[rootExit] -> Value
```

A canonical CFLOBDD must not contain duplicate values in its final value tuple.

If duplicates arise, corresponding root exits must be merged and the merge propagated recursively into the grouping.

---

## 8. Canonical Structural Invariants

The representation is canonical only when all invariants are maintained.

### Invariant 1: A return map is identity and ordered

If the A child has `r` exits, the current grouping has exactly `r` middle vertices and:

```text
AReturnTuple = [1,2,...,r]
```

Thus A exits correspond one-to-one and in order with middle vertices.

Implementation implication: a reduced internal grouping normally does not need an arbitrary A return map. Some repositories still store it for uniformity.

### Invariant 2a: Each B return tuple is injective

Within a single B return tuple, no current-grouping exit may appear twice.

Bad:

```text
[1,1,2]
```

Valid:

```text
[1,3,2]
```

If a B child has two exits that should both map to the same parent exit, the B child itself must first be reduced.

### Invariant 2b: New exits are introduced compactly left-to-right

Across B-connections, exits are numbered in first-occurrence order.

If previous B return tuples have introduced exits `1..M`, then new exit numbers introduced by the next B return tuple must be:

```text
M+1, M+2, ..., M+r
```

in that order.

Existing exits `1..M` may be reused.

Example:

```text
BReturnTuple[1] = [1,2]
BReturnTuple[2] = [2,3]
```

is canonical.

An equivalent numbering such as `[2,1]` followed by `[1,3]` may violate canonical left-to-right numbering.

This invariant explains why the helper is called `CollapseClassesLeftmost`.

### Invariant 3: Equal proto-CFLOBDDs are physically shared

Within a proto-CFLOBDD, there must not be two separate equal lower-level structures.

Enforce via hash-consing/interning.

### Invariant 4: Duplicate `(BConnection, BReturnTuple)` pairs are forbidden

Two middle vertices may point to the same B grouping only if their return tuples differ.

If both grouping and return tuple are equal, the middle vertices are semantically indistinguishable and must be merged. This merge must then propagate backward to the A-connection.

### Invariant 5: Root grouping heads a valid proto-CFLOBDD

The top-level structure must satisfy all grouping invariants recursively.

### Invariant 6: Top-level terminal values are distinct

```text
valueTuple[i] != valueTuple[j] for i != j
```

Duplicate terminal values trigger reduction.

---

## 9. Canonicality, Hash-Consing, and Equality

These are related but distinct concepts.

### 9.1 Canonicality

For a fixed variable ordering, equivalent functions have isomorphic CFLOBDDs.

This is a semantic theorem.

### 9.2 Hash-Consing

Whenever a grouping is constructed:

```cpp
GroupingHandle RepresentativeGrouping(candidate);
```

the unique table is checked:

- return existing equivalent object if present;
- otherwise insert the candidate.

Similarly:

```cpp
CFLOBDDHandle RepresentativeCFLOBDD(grouping, valueTuple);
```

### 9.3 Combined Effect

Canonicality + hash-consing implies:

```cpp
f == g
```

can be implemented as pointer equality.

This is important for:

- fixed-point computations,
- memoization keys,
- identity shortcuts,
- structural comparison,
- reducing repeated recursive work.

### 9.4 Functional/Persistent Architecture

All edges point from level `k` to level `k-1`. The hierarchy is acyclic by level.

Therefore implementation should prefer:

- immutable grouping objects,
- shared handles/smart pointers,
- reference counting or garbage collection,
- no in-place mutation of canonical objects,
- construction followed by interning.

Never mutate an object after inserting it into a unique table.

---

## 10. Caches

Two different hash-based mechanisms are required.

### 10.1 Unique Tables

Purpose:

```text
structurally equal object -> same pointer
```

Typical tables:

- grouping unique table,
- CFLOBDD unique table,
- return tuple unique table,
- value tuple unique table,
- pair tuple / reduction tuple unique tables.

### 10.2 Computed Tables / Function Caches

Purpose:

```text
same operation arguments -> reuse operation result
```

Important caches:

```text
PairProduct(g1,g2)
Reduce(g,reductionTuple)
MatrixMultOnGrouping(g1,g2)
Restriction(...)
ExistentialQuantification(...)
VariableShift(...)
```

Cache keys should use canonical handles where possible.

For commutative operations, consider normalized key ordering only if the returned metadata is symmetric or transformed correctly. `PairProduct(g1,g2)` returns ordered exit pairs, so blindly swapping operands can break semantics.

---

## 11. Interpretation Algorithm

Reference pseudocode:

```cpp
Value InterpretCFLOBDD(const CFLOBDD& n, const Assignment& a) {
    ExitIndex e = InterpretGrouping(n.grouping, a);
    return n.valueTuple[e];
}

ExitIndex InterpretGrouping(
    const GroupingHandle& g,
    const AssignmentSlice& a
) {
    if (g == ForkGrouping) {
        return a[0] ? 1 : 0;
    }

    if (g == DontCareGrouping) {
        return 0;
    }

    if (g == NoDistinctionProtoCFLOBDD(g.level)) {
        return 0;
    }

    size_t half = 1ULL << (g.level - 1);

    auto aA = a.first(half);
    auto aB = a.subspan(half, half);

    ExitIndex middle =
        InterpretGrouping(g.AConnection, aA);

    ExitIndex childExit =
        InterpretGrouping(g.BConnections[middle], aB);

    return g.BReturnTuples[middle][childExit];
}
```

Notes:

- The A result selects a B-connection.
- The B return tuple, not the B child alone, determines the parent exit.
- The no-distinction shortcut avoids traversing an exponentially large conceptual subtree.
- Avoid physically copying assignment halves; use spans/ranges/index intervals.

---

## 12. Primitive Construction

### 12.1 No-Distinction Proto-CFLOBDD

Represents a one-exit partition for all assignments.

```cpp
GroupingHandle NoDistinction(unsigned k) {
    if (k == 0) return DontCareGroupingSingleton;

    InternalGrouping g(k);
    g.AConnection = NoDistinction(k - 1);
    g.AReturnTuple = [0];

    g.BConnections = [g.AConnection];
    g.BReturnTuples = [[0]];
    g.numberOfExits = 1;

    return RepresentativeGrouping(g);
}
```

Only one new grouping is required per level.

### 12.2 Constant Function

```cpp
CFLOBDD Constant(unsigned k, Value v) {
    return RepresentativeCFLOBDD(
        NoDistinction(k),
        [v]
    );
}
```

### 12.3 Projection Function

Construct:

```text
f(x0,...,x_{2^k-1}) = x_i
```

Recurrence:

- If `i` lies in the first half, put the projection recursively in A and use no-distinction B children to preserve the A result.
- If `i` lies in the second half, use no-distinction A and put the projection recursively in the single B child.

Conceptual pseudocode:

```cpp
GroupingHandle ProjectionProto(unsigned k, unsigned i) {
    if (k == 0) return ForkGroupingSingleton;

    unsigned half = 1U << (k - 1);
    InternalGrouping g(k);

    if (i < half) {
        g.AConnection = ProjectionProto(k - 1, i);
        g.AReturnTuple = [0,1];

        g.BConnections = [
            NoDistinction(k - 1),
            NoDistinction(k - 1)
        ];
        g.BReturnTuples = [
            [0],
            [1]
        ];
        g.numberOfExits = 2;
    } else {
        g.AConnection = NoDistinction(k - 1);
        g.AReturnTuple = [0];

        g.BConnections = [
            ProjectionProto(k - 1, i - half)
        ];
        g.BReturnTuples = [
            [0,1]
        ];
        g.numberOfExits = 2;
    }

    return RepresentativeGrouping(g);
}
```

---

## 13. Unary Operations

### 13.1 Value-Tuple-Only Operations

When an operation preserves equality distinctions among terminal values, reuse the grouping and transform only the value tuple.

Examples:

- Boolean complement for a two-terminal diagram,
- scalar multiplication by a nonzero scalar over an integral domain,
- injective relabeling,
- multiplying complex amplitudes by a nonzero phase.

```cpp
return RepresentativeCFLOBDD(
    n.grouping,
    map(n.valueTuple, unaryOp)
);
```

### 13.2 Operations That May Create Duplicate Values

Examples:

```text
square: [-1, 1] -> [1, 1]
absolute value
thresholding
rounding
multiplication by zero
non-injective relabeling
```

Required procedure:

```text
newValues = map(oldValues, op)
(canonicalValues, reductionTuple) =
    CollapseClassesLeftmost(newValues)
newGrouping =
    Reduce(oldGrouping, reductionTuple)
return RepresentativeCFLOBDD(newGrouping, canonicalValues)
```

Do not leave duplicate values in the top-level value tuple.

---

## 14. CollapseClassesLeftmost

This utility canonicalizes equivalence classes by first appearance.

Example:

```text
input:
[2,2,1,1,4,1,1]

projected classes:
[2,1,4]

renumbered classes:
[1,1,2,2,3,2,2]
```

0-based implementation:

```cpp
template<class T>
std::pair<std::vector<T>, std::vector<unsigned>>
CollapseClassesLeftmost(const std::vector<T>& xs) {
    std::unordered_map<T, unsigned> firstClass;
    std::vector<T> projected;
    std::vector<unsigned> reduction;
    reduction.reserve(xs.size());

    for (const T& x : xs) {
        auto [it, inserted] =
            firstClass.emplace(x, projected.size());

        if (inserted) {
            projected.push_back(x);
        }

        reduction.push_back(it->second);
    }

    return {projected, reduction};
}
```

Requirements:

- equality/hash for `T` must be semantically valid;
- floating-point values require an explicit canonical representation;
- iteration order must be deterministic;
- preserve left-to-right first occurrence.

---

## 15. Binary Apply Overview

For two same-level CFLOBDDs `n1`, `n2`:

```text
result(x) = op(n1(x), n2(x))
```

Algorithm:

```cpp
(g, pairTuple) = PairProduct(n1.grouping, n2.grouping);

deducedValues =
    [op(n1.valueTuple[i], n2.valueTuple[j])
     for (i,j) in pairTuple];

(values, reductionTuple) =
    CollapseClassesLeftmost(deducedValues);

gReduced =
    Reduce(g, reductionTuple);

return RepresentativeCFLOBDD(
    gReduced,
    values
);
```

`PairProduct` builds the common refinement of the two operand partitions.

`Reduce` removes distinctions that become unnecessary after `op` is applied.

This separation is fundamental.

---

## 16. PairProduct

### 16.1 Contract

```cpp
PairProduct(g1, g2)
    -> (g, pairTuple)
```

`g` is a proto-CFLOBDD.

Each exit `r` of `g` corresponds to an ordered pair:

```text
pairTuple[r] = (exit of g1, exit of g2)
```

The set of assignments reaching `r` is the intersection of assignments reaching those two operand exits.

### 16.2 Base Cases

#### Both no-distinction

```text
result grouping = either no-distinction operand
pairTuple = [(1,1)]
```

#### Left no-distinction

```text
result grouping = g2
pairTuple = [(1,j) for each exit j of g2]
```

#### Right no-distinction

```text
result grouping = g1
pairTuple = [(i,1) for each exit i of g1]
```

#### Both fork

Only equal bits can occur on the same assignment:

```text
pairTuple = [(false,false), (true,true)]
```

Not the full Cartesian product of four pairs.

### 16.3 Recursive Case

First pair the A-connections:

```cpp
(gA, ptA) =
    PairProduct(g1.AConnection, g2.AConnection);
```

Each pair in `ptA` corresponds to one reachable pair of middle vertices. Therefore create exactly one result B-connection per reachable middle pair.

For each:

```cpp
(i1, i2) = ptA[m];

(gB, ptB) =
    PairProduct(
        g1.BConnections[i1],
        g2.BConnections[i2]
    );
```

For every child exit pair `(j1,j2)` in `ptB`, route through the operand B return tuples:

```cpp
e1 = g1.BReturnTuples[i1][j1];
e2 = g2.BReturnTuples[i2][j2];
```

The ordered pair `(e1,e2)` identifies a result exit.

Maintain `pairTuple` in first-occurrence order. Reuse the existing result-exit index when the pair has already occurred.

### 16.4 Why PairProduct Is Not a Naive Cartesian Product

It constructs only reachable exit pairs induced by the same assignment.

This is analogous to synchronized product construction.

### 16.5 Complexity

With memoization:

- each pair of same-level operand groupings is processed at most once;
- result middle count is bounded by the product of operand middle counts;
- result exit count is bounded by the product of operand exit counts.

Expected cost is polynomial and commonly described as product-bounded in operand sizes.

### 16.6 PairProduct Hazards

- Pair order matters.
- `pairTuple` order affects canonical numbering.
- Do not iterate unordered containers when assigning new exit indices.
- The result passed to `RepresentativeGrouping` must already satisfy grouping invariants 1-4.
- No-distinction shortcuts are important for scalability.
- Cache key must include both grouping handles in order.

---

## 17. Reduce

### 17.1 Contract

```cpp
Reduce(g, reductionTuple) -> g'
```

`reductionTuple` maps every old exit of `g` to a new exit equivalence class.

Example:

```text
old exits:       1 2 3 4 5 6 7
reductionTuple:  1 1 2 2 3 2 2
```

means:

```text
{1,2}     -> new exit 1
{3,4,6,7} -> new exit 2
{5}       -> new exit 3
```

The tuple must already use compact left-to-right class numbering.

### 17.2 Fast Paths

#### Identity reduction

```text
[1,2,...,n]
```

Return `g`.

#### All exits merged

Return `NoDistinctionProtoCFLOBDD(g.level)`.

### 17.3 Core Direction

Reduction propagates backwards:

```text
parent exits
    -> B return tuples
        -> reduce B children
            -> merge duplicate B connection/return pairs
                -> reduce A child
```

This is the opposite conceptual direction from interpretation.

### 17.4 Processing Each B-Connection

For B-connection `i`:

```cpp
deducedReturnClasses =
    [reductionTuple[parentExit]
     for parentExit in g.BReturnTuples[i]];
```

Collapse duplicates:

```cpp
(inducedReturnTuple,
 inducedChildReduction) =
    CollapseClassesLeftmost(deducedReturnClasses);
```

Recursively reduce child:

```cpp
h =
    Reduce(
        g.BConnections[i],
        inducedChildReduction
    );
```

Insert canonical `(h, inducedReturnTuple)` into the new parent:

```cpp
position =
    InsertBConnection(
        newG,
        h,
        inducedReturnTuple
    );
```

Record the position in `reductionTupleA`.

### 17.5 Why `InsertBConnection` Matters

If two reduced B branches become identical in both:

```text
child grouping
return tuple
```

then they must share a single middle vertex.

`InsertBConnection` enforces invariant 4 and returns the existing/new B index.

Example:

```text
old middle 1 -> (H, [1])
old middle 2 -> (H, [1])
```

After reduction both are identical, so new grouping has one middle vertex.

### 17.6 Reducing the A-Connection

After all B connections are processed:

```text
reductionTupleA[old middle]
    = new B-connection index
```

Collapse that tuple:

```cpp
(inducedAReturnTuple,
 inducedAReduction) =
    CollapseClassesLeftmost(reductionTupleA);
```

Then:

```cpp
newA =
    Reduce(
        g.AConnection,
        inducedAReduction
    );

newG.AConnection = newA;
newG.AReturnTuple = inducedAReturnTuple;
```

Finally hash-cons `newG`.

### 17.7 Reduce Pseudocode

```cpp
GroupingHandle Reduce(
    GroupingHandle g,
    ReductionTupleHandle reduction
) {
    if (reduction->isIdentity()) return g;

    if (reduction->numberOfClasses() == 1) {
        return NoDistinction(g->level);
    }

    InternalGroupingBuilder out(g->level);
    out.numberOfExits = reduction->numberOfClasses();

    std::vector<unsigned> reductionA;

    for (size_t i = 0; i < g->BConnections.size(); ++i) {
        std::vector<unsigned> deduced;

        for (unsigned oldParentExit :
             *g->BReturnTuples[i]) {
            deduced.push_back(
                (*reduction)[oldParentExit]
            );
        }

        auto [newReturn, childReduction] =
            CollapseClassesLeftmost(deduced);

        GroupingHandle child =
            Reduce(
                g->BConnections[i],
                InternReductionTuple(childReduction)
            );

        unsigned newMiddle =
            out.InsertBConnection(
                child,
                InternReturnTuple(newReturn)
            );

        reductionA.push_back(newMiddle);
    }

    auto [newAReturn, aReduction] =
        CollapseClassesLeftmost(reductionA);

    out.AConnection =
        Reduce(
            g->AConnection,
            InternReductionTuple(aReduction)
        );

    out.AReturnTuple =
        InternReturnTuple(newAReturn);

    return RepresentativeGrouping(out.Finish());
}
```

### 17.8 Why Reduction Cannot Be Replaced by Hash-Consing Alone

BDD reduction is often achieved locally during bottom-up construction.

CFLOBDD reduction is global across the hierarchy:

- merging root exits can merge exits of B children;
- this can make B branches identical;
- merging B branches changes middle vertices;
- this forces reduction of the A child.

Hash-consing detects equality of already built objects, but does not compute this recursive quotient.

---

## 18. Boolean Operations

Every binary Boolean function can be implemented through Apply because there are only 16 two-input truth tables.

Examples:

```text
AND, OR, XOR, XNOR, implication, NAND, NOR
```

Use a terminal operation object/function:

```cpp
bool op(bool a, bool b);
```

Structural logic remains unchanged.

Useful special cases before generic Apply:

```text
f AND false = false
f AND true  = f
f OR false  = f
f OR true   = true
f XOR false = f
f XOR f     = false
```

Pointer equality makes these shortcuts cheap.

---

## 19. Matrix and Vector Representation

### 19.1 Square Matrix

A `2^n x 2^n` matrix is represented as a function over `2n` bits:

```text
row bits:    x0,...,x_{n-1}
column bits: y0,...,y_{n-1}
```

Common order:

```text
x0, y0, x1, y1, ..., x_{n-1}, y_{n-1}
```

This interleaved order often exposes recursive block structure.

A CFLOBDD level must satisfy:

```text
2^level = 2n
```

so:

```text
level = log2(n) + 1
```

when `n` itself is a power of two in the paper's standard setup.

Repositories may pad dimensions to powers of two.

### 19.2 Vectors

A length-`2^n` vector is represented with `n` Boolean index bits.

Some matrix algorithms convert a vector into a matrix by embedding it as one column and filling the rest with zero.

### 19.3 Variable Order Is Semantically Critical

Check the repository convention before modifying matrix code:

- row-major then column-major,
- interleaved,
- reverse interleaved,
- tensor-level interleaving,
- vocabulary-specific orders such as `VOC12`, `VOC13`, `VOC14`.

Most apparent matrix multiplication bugs in decision-diagram code are actually variable-order or vocabulary-permutation bugs.

---

## 20. Pointwise Matrix Operations vs Matrix Multiplication

Do not confuse:

### Pointwise multiplication

```text
C[i,j] = A[i,j] * B[i,j]
```

This is ordinary binary Apply on terminal values.

### Matrix multiplication

```text
C[i,j] = sum_k A[i,k] * B[k,j]
```

This requires symbolic summation over an index vocabulary and uses a specialized recursive algorithm.

---

## 21. Matrix Multiplication

### 21.1 Main Idea

The grouping structure partitions matrix entries into equivalence classes. Matrix multiplication must compute symbolic weighted sums of pairs of operand exit values.

The algorithm uses a `MatMultTuple`.

Each element is a bilinear polynomial over pairs of operand exits:

```text
bp[(exit1, exit2)] = integer coefficient
```

Interpretation:

```text
sum over (e1,e2):
    coefficient(e1,e2)
    * leftValue[e1]
    * rightValue[e2]
```

### 21.2 Why Coefficients Appear

Multiple internal summation-index assignments can generate the same exit pair. The coefficient counts multiplicity.

### 21.3 Top-Level Algorithm

```text
(g, symbolicTuple) =
    MatrixMultOnGrouping(g1, g2)

concreteValues =
    evaluate each bilinear polynomial using
    n1.valueTuple and n2.valueTuple

(values, reductionTuple) =
    CollapseClassesLeftmost(concreteValues)

g = Reduce(g, reductionTuple)

return RepresentativeCFLOBDD(g, values)
```

### 21.4 Recursive Grouping Algorithm

At each internal level:

1. Recursively symbolically multiply A-connections.
2. The resulting A `MatMultTuple` determines symbolic weighted combinations of B-connection products.
3. For each symbolic A result:
   - recursively multiply required B-connection pairs;
   - transform child bilinear polynomials through the two B return tuples;
   - reduce duplicate symbolic values;
   - create a multi-terminal CFLOBDD whose terminal values are bilinear polynomials;
   - sum weighted symbolic CFLOBDDs using ordinary Apply and scalar multiplication.
4. Collect the resulting symbolic value tuples.
5. Collapse duplicate bilinear polynomials.
6. Reduce the grouping.
7. Return `(grouping, MatMultTuple)`.

### 21.5 Required Algebra

Terminal/value type for symbolic matrix multiplication must support:

- zero,
- addition,
- integer scalar multiplication,
- equality,
- hashing/canonicalization.

Concrete matrix value type must support:

- zero,
- addition,
- multiplication,
- equality,
- hashing or canonical identity.

For floating or complex values, exact equality is problematic. Prefer:

- exact rationals,
- symbolic expressions,
- canonical algebraic numbers,
- controlled quantization with documented semantics,
- high-precision normalized numeric wrappers.

### 21.6 Matrix Multiplication Debugging Checklist

Verify:

- operand levels match;
- variable vocabulary/order matches algorithm assumptions;
- base case matrix size is correct;
- bilinear-polynomial maps combine duplicate supports;
- return tuple lifting transforms both exit indices correctly;
- multiplicities are accumulated, not overwritten;
- symbolic zero entries are canonicalized;
- `CollapseClassesLeftmost` compares full symbolic polynomials;
- final concrete values are reduced;
- cache keys include vocabulary/order mode if relevant.

---

## 22. Path Counting

A grouping can store/cache the number of matched paths reaching each exit.

Base:

```text
DontCareGrouping: [2]
ForkGrouping:     [1,1]
```

Internal recurrence:

```text
count[parentExit] +=
    countA[middle]
    * countB_i[childExit]
```

routed through `BReturnTuple[i][childExit]`.

Counts become enormous because a level-`k` grouping conceptually covers `2^(2^k)` assignments. Use arbitrary-precision integers.

Applications:

- count satisfying assignments,
- normalize weighted distributions,
- sample assignments,
- compute quantum measurement probabilities from squared amplitudes,
- sanity-check partition coverage.

---

## 23. Sampling

For nonnegative terminal weights:

```text
probability(assignment)
  proportional to
terminalWeight(assignment)
```

Sampling requires:

1. path counts per exit,
2. total weighted mass,
3. recursively choosing A and B branches with probabilities proportional to downstream mass.

For complex quantum amplitudes, convert to nonnegative weights such as squared magnitudes before probabilistic sampling.

---

## 24. Complexity Intuition

Let `|G|` denote stored CFLOBDD size, including grouping structures and tuples.

Typical expectations:

- interpretation: `O(2^level)` decisions without no-distinction shortcuts; often much less structurally;
- constant/projection construction: `O(level)`;
- value-tuple relabeling: `O(number of root exits)`;
- PairProduct: polynomial, commonly bounded by product of operand structural sizes;
- Reduce: polynomial in grouping and reduction tuple size, with memoization;
- Binary Apply: product construction plus reduction;
- equality: expected `O(1)` pointer comparison under canonical hash-consing;
- matrix multiplication: substantially more expensive and strongly dependent on symbolic tuple growth.

Do not confuse represented variable count with stored size:

```text
variables = 2^level
assignment space = 2^(2^level)
```

A small level increase doubles variables and squares assignment-space size.

---

## 25. Repository Architecture Guidance

A maintainable CFLOBDD repository should have clear layers.

### Layer 1: Canonical immutable data

- grouping classes/handles,
- return tuple handles,
- value tuple handles,
- pair tuples,
- reduction tuples,
- bilinear polynomial objects.

### Layer 2: Unique tables and memory management

- grouping interning,
- CFLOBDD interning,
- tuple interning,
- reference counting / smart pointers,
- unique-table cleanup.

### Layer 3: Structural algorithms

- no-distinction,
- projection,
- PairProduct,
- Reduce,
- InsertBConnection,
- vocabulary transforms,
- restriction/quantification.

### Layer 4: Generic value-domain operations

- terminal Apply,
- scalar multiplication,
- unary map,
- arithmetic wrappers.

### Layer 5: Matrix/quantum algorithms

- identity,
- Kronecker product,
- matrix multiplication,
- vector conversion,
- gate constructors,
- state evolution,
- measurement/sampling.

### Layer 6: Tests and diagnostics

- invariant checker,
- semantic evaluator for small levels,
- node/edge/exit statistics,
- cache statistics,
- DOT printer,
- matrix/vector dump for small dimensions.

Avoid mixing terminal arithmetic with grouping canonicalization.

---

## 26. Invariant Checker

Before modifying algorithms, implement or locate a recursive validator.

Suggested checks:

```cpp
ValidateGrouping(g):
    assert(g.level >= 0)

    if level == 0:
        assert(g is canonical fork or don't-care singleton)
        return

    assert(AConnection.level == g.level - 1)
    assert(BConnections.size == AConnection.numberOfExits)
    assert(BReturnTuples.size == BConnections.size)
    assert(AReturnTuple == identity)
    assert(g.numberOfExits > 0)

    for each B i:
        assert(BConnections[i].level == g.level - 1)
        assert(BReturnTuples[i].size ==
               BConnections[i].numberOfExits)
        assert(BReturnTuples[i] is injective)
        assert(all targets < g.numberOfExits)

    assert exits are introduced compactly in
           first-occurrence order)

    assert no duplicate
           (BConnection, BReturnTuple) pairs)

    recursively validate children
```

For a full CFLOBDD:

```cpp
assert(valueTuple.size ==
       grouping.numberOfExits)

assert(valueTuple has no duplicates)
```

During intermediate algorithms, distinguish:

- mock grouping,
- valid proto-CFLOBDD,
- reduced/canonical full CFLOBDD.

Do not apply full-result invariants prematurely to intermediate structures.

---

## 27. Testing Strategy

### 27.1 Small-Level Exhaustive Semantic Tests

For levels 0-3:

- enumerate every assignment;
- compare CFLOBDD interpretation against a truth table or dense matrix;
- test all binary Boolean operations;
- test Reduce against explicit quotienting;
- test PairProduct exit-pair semantics.

### 27.2 Canonicality Tests

Construct the same function via different expressions:

```text
x AND y
y AND x
NOT(NOT x)
x OR false
x XOR x
```

Expected:

```text
semantic equality -> same canonical CFLOBDD pointer
```

when variable ordering and terminal representation match.

### 27.3 Reduction Tests

Cases:

- identity reduction,
- all exits merge,
- only one B child reduces,
- two B connections become identical,
- A middle vertices merge,
- nested propagation across multiple levels.

### 27.4 Matrix Tests

For small dimensions compare with dense arithmetic:

- identity multiplication,
- zero matrix,
- permutation matrices,
- Hadamard-like matrices,
- associativity within exact arithmetic,
- Kronecker identities,
- vector-matrix multiplication.

### 27.5 Hash-Consing Tests

- repeated construction returns identical handles;
- unique-table size stabilizes;
- no canonical object mutates;
- reference-count cleanup does not leave dangling cache entries.

---

## 28. Common Bugs

### 28.1 1-Based vs 0-Based Indexing

Paper pseudocode uses 1-based indices. Most C++ uses 0-based.

Audit:

- exit indices,
- middle indices,
- return tuples,
- projection index,
- pair tuples,
- reduction tuple class IDs.

### 28.2 Duplicate Terminal Values

Returning:

```text
grouping + [v,v]
```

as a final CFLOBDD violates canonicality.

Must collapse and reduce.

### 28.3 Non-Injective B Return Tuples

A parent B return tuple cannot merge two child exits directly. Reduce the child first.

### 28.4 Duplicate B Connections

Equal child pointer plus equal return tuple means duplicate middle vertices. Use `InsertBConnection`.

### 28.5 Unstable Ordering

Using unordered iteration to assign exit numbers breaks deterministic canonical numbering.

### 28.6 Mutating Interned Objects

Never modify canonical grouping/tuple objects after interning.

### 28.7 Cache Lifetime

A computed-table entry may retain or refer to objects whose unique-table lifecycle has ended. Make ownership policy explicit.

### 28.8 Floating-Point Equality

Canonicalization requires an equivalence relation compatible with hashing. Raw approximate floating-point equality is unsafe.

### 28.9 Incorrect No-Distinction Shortcut

Only use no-distinction when all assignments truly reach one exit. Do not confuse a one-terminal value tuple with a one-exit grouping before reduction.

### 28.10 Matrix Vocabulary Mismatch

Pointwise-correct code can still produce a wrong matrix if row/column/summation bits are misordered.

---

## 29. How an AI Agent Should Explore a CFLOBDD Repository

Recommended order:

1. Locate core handle/node/grouping classes.
2. Identify unique tables and representative constructors.
3. Find return tuple and value tuple types.
4. Determine index convention.
5. Locate no-distinction and fork singletons.
6. Read the interpretation/evaluation function.
7. Read invariant or debug-printing code.
8. Read PairProduct.
9. Read Reduce and InsertBConnection together.
10. Read generic Apply wrapper.
11. Read matrix vocabulary/order definitions.
12. Read matrix multiplication only after the above.
13. Run small exhaustive tests before editing optimization code.

Useful repository searches:

```bash
rg -n "PairProduct|CrossProduct"
rg -n "Reduce\\("
rg -n "InsertBConnection"
rg -n "NoDistinction|DontCare|Fork"
rg -n "ReturnTuple|ReturnMap"
rg -n "Representative|Canonical|UniqueTable|HashCons"
rg -n "valueTuple|ValuesList"
rg -n "MatrixMultiply|MatMult|Bilinear"
rg -n "VOC1|VOC2|interleav|vocabulary"
```

Before making a change, the agent should state:

```text
- which invariant is affected;
- whether output is mock, proto, or full CFLOBDD;
- whether canonical numbering is preserved;
- whether a unique-table insertion is safe;
- whether a computed-table key remains valid;
- whether terminal-value equality is exact.
```

---

## 30. Modification Protocol

For any structural change:

### Step 1: State semantic contract

Example:

```text
Reduce(g,r) returns a grouping whose assignments reach
new exit r[oldExit], preserving the represented partition quotient.
```

### Step 2: State preserved invariants

At minimum:

```text
A identity return map
B return injectivity
compact left-to-right numbering
no duplicate proto-CFLOBDDs
no duplicate (B child, return tuple)
```

### Step 3: Identify intermediate form

Is the object:

- arbitrary/mock,
- invariant-valid proto,
- fully reduced grouping,
- full CFLOBDD with distinct terminal values?

### Step 4: Add exhaustive tests at low levels

Do not rely only on large quantum benchmarks.

### Step 5: Measure structural effects

Report:

- groupings by level,
- exits/middles,
- unique-table entries,
- cache hits/misses,
- peak memory,
- terminal count.

A runtime speedup that destroys sharing may be a regression at larger levels.

---

## 31. Compact Algorithm Reference

### Interpret

```text
A on first half
B selected by A exit on second half
route B exit through B return tuple
```

### Constant

```text
NoDistinction(level) + [value]
```

### Projection

```text
recurse into A if variable is in first half
otherwise recurse into B
```

### Unary map

```text
map values
collapse duplicates
Reduce if necessary
```

### Binary Apply

```text
PairProduct
apply op to exit-value pairs
CollapseClassesLeftmost
Reduce
RepresentativeCFLOBDD
```

### PairProduct

```text
synchronized common refinement of assignment partitions
returns result grouping + ordered operand-exit pairs
```

### Reduce

```text
merge root exits
propagate merge into B children
merge duplicate reduced B branches
propagate middle merge into A child
```

### Matrix multiplication

```text
recursive symbolic products
bilinear-polynomial terminal descriptors
evaluate descriptors at root
collapse and reduce
```

---

## 32. Key Design Rules

1. A grouping represents a partition of assignments, not terminal values.
2. The value tuple labels root partition classes.
3. Variables are determined by traversal context, not leaf identity.
4. Legal paths obey call/return matching.
5. A level-`k` grouping consumes `2^k` variables.
6. A exits correspond identically to middle vertices.
7. B return tuples are injective.
8. Exit numbering is first-occurrence, left-to-right.
9. Equal lower-level structures must be shared.
10. Duplicate `(B child, return tuple)` pairs imply mergeable middle vertices.
11. Duplicate root values imply mergeable root exits.
12. PairProduct creates distinctions; Reduce removes unnecessary distinctions.
13. Hash-consing is not a substitute for Reduce.
14. Canonical objects are immutable.
15. Deterministic ordering is part of correctness.
16. Exact terminal equality is part of canonicality.
17. Matrix variable order is part of semantics.
18. Test structural invariants and dense semantics independently.

---

## 33. Final Perspective

The key conceptual transition is:

```text
BDD:
    variable-labeled DAG with suffix sharing

CFLOBDD:
    hierarchical call/return program whose exits encode
    equivalence classes of assignments
```

The most important implementation transition is:

```text
BDD Apply:
    recursively combine nodes, reduce locally

CFLOBDD Apply:
    recursively construct a common partition,
    compute terminal equivalences,
    then propagate the quotient backwards through
    B-connections and A-connections
```

An AI agent that understands only the public `CFLOBDD` wrapper but not:

- return tuples,
- left-to-right exit numbering,
- PairProduct metadata,
- Reduce propagation,
- unique-table discipline,

is likely to produce code that passes small examples but violates canonicality or loses compression.

Treat the structural invariants as part of the type system of the implementation, even if the programming language does not encode them statically.
