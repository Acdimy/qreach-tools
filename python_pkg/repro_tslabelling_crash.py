#!/usr/bin/env python3
r"""Minimal, stable reproduction of SIGSEGV/SIGBUS in tsLabelling.

Reproduces: pe_7 with Pauli-Y injection → crash at ``tsLabelling()``.

Usage::

    cd python_pkg && ../.venv/bin/python repro_tslabelling_crash.py

Expected: exit code -10 (SIGBUS) or -11 (SIGSEGV) at ``tsLabelling``.

Why a subprocess?  CFLOBDD global state survives across calls to
``initializeTransitionSystem()``, so the crash is only reproducible in a
**fresh** OS process.  This script spawns itself in a child process to
guarantee that.

The exact parameters are derived from the VeriQBench experiment
(seed=1, pe_7 at file_index=5):

* QASM:        ``benchmark/pe/pe_7.qasm`` (8 qubits, 42 gates)
* Initial:     ``01010001``  (deterministic from seed=3010007)
* Injection:   ``Y`` gate at instruction position 36 on qubit 1
"""

import subprocess
import sys
from pathlib import Path

_RUNNER = r'''
import copy, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import pyqreach
from qiskit import QuantumCircuit
from qiskit.circuit import CircuitInstruction
from qiskit.circuit.library import YGate
from parse_qiskit import parse_qiskit_cir_lazy
from qctl import tsLabelling

pyqreach.initializeTransitionSystem()

qc_orig = QuantumCircuit.from_qasm_file("benchmark/pe/pe_7.qasm")

# Inject Y at pos 36 on qubit 1
qc_inj = copy.deepcopy(qc_orig)
qc_inj.data.insert(36, CircuitInstruction(YGate(), [qc_inj.qubits[1]], []))

print("[1] Parse+fp injected ...", flush=True)
ts = pyqreach.TransitionSystem()
r = parse_qiskit_cir_lazy(qc_inj, 8, ts, initial_state="01010001", return_metadata=True)
ts.computingFixedPointPost()

print("[2] Parse+fp clean ...", flush=True)
ts2 = pyqreach.TransitionSystem()
r2 = parse_qiskit_cir_lazy(qc_orig, 8, ts2, initial_state="01010001", return_metadata=True)
ts2.computingFixedPointPost()
expected = ts2.Locations[r2.result_locations[0]].lowerBound

print("[3] tsLabelling (CRASH SITE) ...", flush=True)
for loc in r.result_locations:
    ts.setLabel(loc, "final")
tsLabelling(ts, expected, "debug_op", locList=list(r.result_locations))
print("NO CRASH (unexpected)", flush=True)
'''


def main() -> None:
    me = Path(__file__).resolve()
    code = _RUNNER.replace(
        'str(Path(__file__).resolve().parent)',
        repr(str(me.parent)),
    )

    print("=== tsLabelling SIGSEGV repro  (pe_7 + Y injection) ===", flush=True)
    proc = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(me.parent),
        timeout=60,
    )
    rc = proc.returncode
    if rc == -11 or rc == 139 or rc == 245:
        print(f"CRASH: SIGSEGV (exit {rc})  — reproduced successfully", flush=True)
    elif rc == -10 or rc == 138:
        print(f"CRASH: SIGBUS  (exit {rc})  — reproduced successfully", flush=True)
    elif rc == 0:
        print("WARNING: no crash this time (Heisenbug — try again)", flush=True)
    else:
        print(f"UNEXPECTED exit {rc}", flush=True)


if __name__ == "__main__":
    main()
