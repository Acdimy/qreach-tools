#!/usr/bin/env python3
"""
Backend cross-check: QReach post-image vs Qiskit `Statevector` oracle.

This is the "dense cross-check oracle" called for by the backend contract
(`docs/agent-handoffs/backend-replacement-api-contract.md`, Phase 6 item 2):
for a matrix of small *unitary* circuits, verify that QReach's fixed-point
post image equals Qiskit's `Statevector` evolution up to global phase. It is
the safety net for the "accept double precision, rely on a Qiskit oracle"
decision, and it catches variable-order / endianness mismatches the moment a
new DD backend (CFLOBDD or LimTDD) is dropped in.

Backend-agnostic: it runs identically under CFLOBDD or LimTDD (only the
`pyqreach` build differs).

Scope: unitary circuits (deterministic pure states) only. Circuits with
measurements / resets / control flow produce multiple leaves or mixed states
and are covered elsewhere (`test_RUS.py`, `test_lazy_measurement.py`).

Run from `python_pkg/`:
    ../.venv/bin/python workflow_tests/test_backend_crosscheck.py
(or) pytest workflow_tests/test_backend_crosscheck.py
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# IMPORTANT: import Qiskit BEFORE pyqreach to avoid CFLOBDD conflicts.
from qiskit import QuantumCircuit  # noqa: E402
from qiskit.quantum_info import Statevector  # noqa: E402
import numpy as np  # noqa: E402

import pyqreach  # noqa: E402
from qreach import symbolic_available  # noqa: E402
from qreach.parse_qiskit import parse_qiskit_cir  # noqa: E402
from qreach.qctl import set_initial_state  # noqa: E402

# Known backend bugs the cross-check exercises. Keyed by (circuit name -> reason).
# These fail under LimTDD because `DDMatrix::MatrixMultiply` (matrix x matrix)
# is not yet compactified / correct (see
# docs/agent-handoffs/limtdd-matrixmultiply-bug-report.md). Under CFLOBDD
# (MatrixMultiplyV4 is mature) they should pass.
KNOWN_LIMTDD_FAILURES = {
    # CSX = H(target)·CP(π/2)·H(target) — needs MatrixMultiply even for ctrl<tgt.
    "csx": "LimTDD MatrixMultiply(matrix×matrix) throws 'size mismatch' (H·C·H composition)",
    # CX with ctrl>tgt goes through S·C·S SWAP conjugation via MatrixMultiply.
    "nonadj_cx": "LimTDD MatrixMultiply silently mis-computes the S·C·S SWAP conjugation (CX ctrl>tgt)",
}


# ---------------------------------------------------------------------------
# Oracle: Qiskit Statevector -> QReach QOperation
# ---------------------------------------------------------------------------

def qiskit_amplitudes(qc: QuantumCircuit, num_qubits: int, init_qreach: str):
    """Evolve |init> under `qc` with Qiskit, return split (re, im) amplitude lists.

    Pure Qiskit — constructs no `pyqreach` object, so it is safe to call before
    the backend is initialized.

    Conventions (verified against the existing simulation tests):
      - Qiskit label is little-endian: leftmost char = qubit (n-1), rightmost = qubit 0.
      - QReach amplitude index is big-endian: qubit 0 = MSB.
      So a Qiskit label maps to a QReach basis index via string reversal.
      - QReach `QOperation(amps, qubits)` takes SPLIT format:
        [re_0 .. re_{N-1}, im_0 .. im_{N-1}], N = 2^num_qubits.
    """
    init_qiskit = init_qreach[::-1]
    sv = Statevector.from_label(init_qiskit).evolve(qc)
    n_states = 1 << num_qubits
    re = [0.0] * n_states
    im = [0.0] * n_states
    for label, amp in sv.to_dict().items():
        basis = int(str(label)[::-1], 2)
        re[basis] = float(np.real(amp))
        im[basis] = float(np.imag(amp))
    return re, im


def check_circuit(name, num_qubits, build, init_qreach):
    """Run one circuit through Qiskit (oracle) and QReach, assert they agree."""
    qc = QuantumCircuit(num_qubits)
    build(qc)

    # Compute the Qiskit oracle amplitudes up front (pure Qiskit, no pyqreach).
    re, im = qiskit_amplitudes(qc, num_qubits, init_qreach)

    # Create the TransitionSystem FIRST: its constructor initializes the DD
    # backend. Constructing a `QOperation` from raw amplitudes before this point
    # corrupts the (LimTDD) backend and makes post-image propagation yield an
    # empty (lower_dim=0) leaf.
    ts = pyqreach.TransitionSystem()
    parse_qiskit_cir(qc, num_qubits, ts)
    expected = pyqreach.QOperation(re + im, num_qubits)
    set_initial_state(ts, init_qreach)
    ts.computingFixedPointPost()

    leaf_locs = [i for i in range(ts.getLocationNum()) if ts.isLeafLoc(i)]
    assert leaf_locs, f"{name}: no leaf locations"

    for loc in leaf_locs:
        upper_dim, lower_dim = ts.printDims(loc)
        assert lower_dim == 1, (
            f"{name}: leaf {loc} should be a pure state, got lower_dim={lower_dim}"
        )
        assert ts.satisfy(loc, expected), (
            f"{name}: leaf {loc} state != Qiskit Statevector (dims=({upper_dim},{lower_dim}))"
        )
    print(f"  {name:<28} OK  ({len(leaf_locs)} leaf, {ts.getLocationNum()} locs)")
    return True


# ---------------------------------------------------------------------------
# Circuit matrix: one builder per gate type / ordering / entanglement shape.
# ---------------------------------------------------------------------------

# (name, num_qubits, build(qc), init_qreach)  — init in QReach convention.
CIRCUITS = [
    ("h_on_q0", 1, lambda qc: qc.h(0), "0"),
    ("x_on_q0", 1, lambda qc: qc.x(0), "0"),
    ("y_on_q0", 1, lambda qc: qc.y(0), "0"),
    ("z_on_q1", 2, lambda qc: qc.z(1), "0" * 2),
    ("s_t_on_q0", 1, lambda qc: (qc.s(0), qc.t(0)), "0"),
    ("sdg_tdg", 1, lambda qc: (qc.sdg(0), qc.tdg(0)), "0"),
    ("bell_h_cx01", 2, lambda qc: (qc.h(0), qc.cx(0, 1)), "00"),
    ("bell_h_cx10", 2, lambda qc: (qc.h(0), qc.cx(1, 0)), "00"),
    ("ghz3", 3, lambda qc: (qc.h(0), qc.cx(0, 1), qc.cx(1, 2)), "000"),
    ("toffoli", 3, lambda qc: (qc.x(0), qc.x(1), qc.ccx(0, 1, 2)), "000"),
    ("swap", 2, lambda qc: (qc.x(0), qc.swap(0, 1)), "00"),
    ("cz", 2, lambda qc: (qc.h(0), qc.h(1), qc.cz(0, 1)), "00"),
    ("cp_angle", 2, lambda qc: (qc.h(0), qc.h(1), qc.cp(0.37, 0, 1)), "00"),
    ("csx", 2, lambda qc: (qc.h(0), qc.csx(0, 1)), "00"),
    ("iswap", 2, lambda qc: (qc.h(0), qc.iswap(0, 1)), "00"),
    ("u3_float", 1, lambda qc: qc.u(0.6, 0.4, 0.9, 0), "0"),
    ("u3_two_qubit", 2, lambda qc: (qc.u(0.6, 0.4, 0.9, 0), qc.u(1.1, -0.7, 0.2, 1)), "00"),
    ("ry_float", 1, lambda qc: qc.ry(0.93346815, 0), "0"),
    ("rz_phase", 1, lambda qc: qc.rz(1.234, 0), "0"),  # global-phase-tolerant
    ("reorder_3q", 3, lambda qc: (qc.h(0), qc.cx(0, 2), qc.swap(0, 1), qc.cx(1, 2), qc.t(0)), "000"),
    ("nonadj_cx", 3, lambda qc: (qc.h(0), qc.cx(0, 2), qc.cx(2, 1)), "000"),
    ("init_not_zero", 2, lambda qc: (qc.h(0), qc.cx(0, 1)), "10"),
    ("init_superpos", 2, lambda qc: (qc.cx(0, 1), qc.z(1)), "+0"),
]


def test_backend_crosscheck():
    is_limtdd = not symbolic_available()
    failures = []
    xfails = []
    for name, n, build, init in CIRCUITS:
        try:
            check_circuit(name, n, build, init)
        except Exception as exc:  # noqa: BLE001 — report and continue
            if is_limtdd and name in KNOWN_LIMTDD_FAILURES:
                xfails.append(name)
                print(f"  {name:<28} XFAIL ({KNOWN_LIMTDD_FAILURES[name]})")
            else:
                failures.append((name, exc))
                print(f"  {name:<28} FAIL: {exc}")
    if failures:
        raise AssertionError(f"{len(failures)}/{len(CIRCUITS)} circuits failed: "
                             f"{[f[0] for f in failures]}")
    print(f"\nAll {len(CIRCUITS)} cross-check circuits PASSED"
          + (f" ({len(xfails)} known-LimTDD XFAIL: {xfails})" if xfails else ""))


if __name__ == "__main__":
    test_backend_crosscheck()
