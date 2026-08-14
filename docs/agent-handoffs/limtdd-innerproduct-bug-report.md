# LimTDD `InnerProduct` scaling bug — breaks Gram-Schmidt / disjunction

**For:** LimTDD backend agent
**From:** QReach side
**Date:** 2026-08-14
**Severity:** correctness — silently wrong results (no crash, no assertion)

## Summary

`DDVector::InnerProduct(a, b)` returns **wrong scalar values** for vectors of
≥ 2 qubits. The error is a scale of `2^(n-1)` (n = qubit count), plus an extra
non-uniform factor on superposition terms. This corrupts the projection
coefficients that QReach's `disjunction` / `span_qops` (Gram-Schmidt
orthogonalization) rely on, so **linearly-dependent vectors are not collapsed**.
The visible symptom is extra support vectors in the QReach RUS workflow
(while_loop + reset + measure), but the root cause is entirely inside the
LimTDD `InnerProduct` primitive — no QReach environment is needed to see it.

The bug is **not** a threshold/precision issue: the values are off by integer
powers of two, e.g. a unit vector's self-inner-product is `8` instead of `1`.

## Minimal reproduction (LimTDD-only, no QReach)

Drop this into `test/` and link against the three core objects. It uses only
`dd/backend/DDVector.hpp`.

```cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dd/backend/DDVector.hpp"

int main() {
    DDVector::Initialize();

    auto ip = [](const char* name, DDVector::DDComplex c, double re, double im) {
        std::printf("  %-22s got (% .6f, % .6f)   expect (% .6f, % .6f)%s\n",
                    name, c.real(), c.imag(), re, im,
                    (std::abs(c.real() - re) < 1e-6 && std::abs(c.imag() - im) < 1e-6)
                        ? "   OK" : "   <<< MISMATCH");
    };

    // <basis|basis> must ALWAYS be 1.
    ip("<|0>|0> (1q)",     DDVector::InnerProduct(DDVector::MkBasisVector(0, 0),
                                                  DDVector::MkBasisVector(0, 0)), 1.0, 0.0);
    ip("<|00>|00> (2q)",   DDVector::InnerProduct(DDVector::MkBasisVector(1, 0),
                                                  DDVector::MkBasisVector(1, 0)), 1.0, 0.0);
    ip("<|1010>|1010> (4q)", DDVector::InnerProduct(DDVector::MkBasisVector(2, 10),
                                                    DDVector::MkBasisVector(2, 10)), 1.0, 0.0);

    const double s2  = 1.0 / std::sqrt(2.0);
    const double s10 = 1.0 / std::sqrt(10.0);
    auto mk = [](const std::vector<double>& a) { return DDVector::InitializeWithAmplitudes(4, a); };
    std::vector<double> a1(32, 0.0), a2(32, 0.0), a3(32, 0.0);
    a1[10] = s2; a1[16 + 10] = s2;                        // v1 = e^{i pi/4} |1010>
    a2[16 + 10] = s2; a2[16 + 14] = -s2;                  // v2 = (i/sqrt2)(|1010>-|1110>)
    a3[16 + 10] = s10; a3[16 + 14] = -3.0 * s10;          // v3 = (i/sqrt10)|1010>-(3i/sqrt10)|1110>
    auto v1 = mk(a1), v2 = mk(a2), v3 = mk(a3);

    ip("<v1|v2>", DDVector::InnerProduct(v1, v2), 0.5, 0.5);
    ip("<v1|v3>", DDVector::InnerProduct(v1, v3), 0.22360679774997896, 0.22360679774997896);
    ip("<v2|v3>", DDVector::InnerProduct(v2, v3), 2.0 / std::sqrt(5.0), 0.0);
    return 0;
}
```

### Actual output

```
== <basis|basis> should ALWAYS be 1 ==
  <|0>|0> (1 qubit)      got ( 1.000000,  0.000000)   expect ( 1.000000,  0.000000)   OK
  <|00>|00> (2 qubit)    got ( 2.000000,  0.000000)   expect ( 1.000000,  0.000000)   <<< MISMATCH
  <|1010>|1010> (4qb)    got ( 8.000000,  0.000000)   expect ( 1.000000,  0.000000)   <<< MISMATCH

== projection coefficients (Gram-Schmidt depends on these) ==
  <v1|v2>                got ( 4.000000,  4.000000)   expect ( 0.500000,  0.500000)   <<< MISMATCH
  <v1|v3>                got (-5.366563, -5.366563)   expect ( 0.223607,  0.223607)   <<< MISMATCH
  <v2|v3>                got ( 7.155418, -0.000000)   expect ( 0.894427,  0.000000)   <<< MISMATCH
```

## Scaling law

`<basis|basis>` (a unit vector with itself, which must be exactly 1):

| qubits | got | expected | factor |
|---|---|---|---|
| 1 | 1 | 1 | 2⁰ ✓ |
| 2 | 2 | 1 | 2¹ |
| 4 | 8 | 1 | 2³ |

So the self-inner-product is scaled by `2^(n-1)` for n qubits. The
superposition cases show a further non-uniform factor: `<v1|v3>` is off by
`-24` (= `-3 · 2³`) while `<v1|v2>` and `<v2|v3>` are off by exactly `2³`.
Because the scale is **not uniform across pairs**, the Gram-Schmidt ratio
`⟨vᵢ|v⟩ / ⟨vᵢ|vᵢ⟩` comes out wrong and the residual is never driven to zero.

## Why this breaks disjunction (QReach symptom, for context)

QReach's `disjunction` / `span_qops` orthogonalizes a new vector against the
existing orthonormal basis using `InnerProduct` coefficients, then drops the
vector if the residual is (approximately) zero. With wrong coefficients:

- **Parallel (collinear) vectors collapse fine** — a single projection onto one
  basis vector, and the ratio happens to still work out.
- **A vector that is a non-trivial linear combination of two basis vectors does
  NOT collapse** — the two projection coefficients are both wrong, the residual
  is non-zero, and the vector is kept.

Concretely, `span_qops([v1, v2, v3])` with `v3 = (-0.4472-0.4472i)·v1 +
1.3416·v2` (linearly dependent) returns **3** vectors on LimTDD but **2** on the
CFLOBDD backend (correct). The extra vector then propagates through the RUS
fixed-point, giving 4 support vectors where CFLOBDD has 2.

## Why the existing `test_ddvector.cpp` misses it

The current `InnerProduct` assertions only cover:

- 1-qubit cases (`<0|0>`, `<0|1>`, `<+|+>`, `<+|->`, `<iphase|iphase>`), where
  the scale factor `2^(n-1)` is `1`, and
- a 2-qubit **orthogonal** case (`<00+11 | 00-11> = 0`), where the answer is
  zero regardless of scaling.

There is no assertion that `2+`-qubit `<v|v> == 1`, which is where the bug
first appears. Adding the `<|00>|00> == 1` check above is the smallest
regression guard.

## Root-cause hypothesis (for your investigation)

`InnerProduct` is implemented as:

```cpp
conj(a) rebuilt densely via GetNonZeroAmplitudes -> InitializeWithAmplitudes,
then Package::cont(conjA, b), then ExtractSingleAmplitude (= detail::sumRemaining).
```

The `2^(n-1)` scale is consistent with `sumRemaining` summing over residual
nodes that `cont` leaves behind (the existing comment in `sumRemaining` already
notes that `cont` does not always collapse to a single terminal). If `cont`
leaves `2^(n-1)` branches for an n-qubit scalar result, summing them all gives
`2^(n-1)` × the correct value. The extra `-3` factor on `<v1|v3>` suggests the
residual structure also depends on the specific non-zero entries of the
operands (not just qubit count), so the correction must make `cont`/`sumRemaining`
collapse to a true scalar, not just divide by a fixed power of two.

Please fix `InnerProduct` (and add the `2+`-qubit `<v|v> == 1` regression test).
After that, the QReach RUS workflow should produce the same 2-vector lower bound
as the CFLOBDD backend.
