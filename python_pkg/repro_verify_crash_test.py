#!/usr/bin/env python3
"""Test which of the 9 crash cases crash at tsLabelling with verify enabled.

Runs each case in its own spawned process (to survive SIGSEGV).
"""
import multiprocessing as mp
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

CASES = [
    ("pe_7",  7,  5,  8),
    ("pe_8",  8,  6,  9),
    ("pe_9",  9,  7,  10),
    ("pe_10", 10, 8,  11),
    ("pe_11", 11, 9,  12),
    ("pe_12", 12, 10, 13),
    ("qft_3", 3,  1,  3),
    ("qft_7", 7,  5,  7),
    ("qft_8", 8,  6,  8),
]


def worker(stem, n, file_index, num_qubits, family):
    """Run full pipeline WITH verify. Dies with SIGSEGV on crash."""
    import copy, random
    from time import perf_counter

    import pyqreach
    from qiskit import QuantumCircuit
    from qiskit.circuit import CircuitInstruction
    from qiskit.circuit.library import XGate, YGate, ZGate
    from parse_qiskit import parse_qiskit_cir_lazy
    from qctl import modelChecking, span_qops, tsLabelling

    SEED = 1
    PY_PKG = Path(__file__).resolve().parent if __file__ else Path('.')
    QASM = PY_PKG / "benchmark" / family / f"{stem}.qasm"
    NUSMV = str((PY_PKG.parent.parent / "NuSMV-2.7.0-macos-universal" / "bin" / "NuSMV").resolve())

    # Init state
    import re
    if family == "pe":
        rng_seed = SEED * 10000 + n + 3_000_000
    else:
        rng_seed = SEED * 10000 + n + 2_000_000
    init_state = "".join(random.Random(rng_seed).choice(["0", "1"]) for _ in range(num_qubits))

    # Load + inject
    qc_orig = QuantumCircuit.from_qasm_file(str(QASM))
    rng = random.Random(SEED + file_index)
    qc_inj = copy.deepcopy(qc_orig)
    pos = rng.randint(0, len(qc_inj.data))
    q = rng.randint(0, qc_inj.num_qubits - 1)
    gate_name = rng.choice(["x", "y", "z"])
    gate_map = {"x": XGate(), "y": YGate(), "z": ZGate()}
    qc_inj.data.insert(pos, CircuitInstruction(gate_map[gate_name], [qc_inj.qubits[q]], []))

    print(f"  [{stem}] init={init_state}  inject={gate_name}@pos{pos}/q{q}", flush=True)

    # Step 1: parse + fp injected
    print(f"  [{stem}] [1] parse+fp injected ...", flush=True)
    ts = pyqreach.TransitionSystem()
    r1 = parse_qiskit_cir_lazy(qc_inj, qc_inj.num_qubits, ts,
                                initial_state=init_state, return_metadata=True)
    ts.computingFixedPointPost()
    print(f"  [{stem}] [1] DONE ({ts.getLocationNum()} locs)", flush=True)

    # Step 2: parse + fp clean
    print(f"  [{stem}] [2] parse+fp clean ...", flush=True)
    ts2 = pyqreach.TransitionSystem()
    r2 = parse_qiskit_cir_lazy(qc_orig, qc_orig.num_qubits, ts2,
                                initial_state=init_state, return_metadata=True)
    ts2.computingFixedPointPost()
    if len(r2.result_locations) == 1:
        expected = ts2.Locations[r2.result_locations[0]].lowerBound
    else:
        expected = span_qops([ts2.Locations[l].lowerBound for l in r2.result_locations])
    print(f"  [{stem}] [2] DONE", flush=True)

    # Step 3: tsLabelling  <-- CRASH SITE
    print(f"  [{stem}] [3] tsLabelling ...", flush=True)
    for loc in r1.result_locations:
        ts.setLabel(loc, "final")
    tsLabelling(ts, expected, "debug_op", locList=list(r1.result_locations))
    print(f"  [{stem}] [3] DONE", flush=True)

    # Step 4: NuSMV
    print(f"  [{stem}] [4] NuSMV ...", flush=True)
    mc = modelChecking(ts, "AG (debug_op <-> final)", nusmv_path=NUSMV)
    print(f"  [{stem}] [4] DONE  sat={mc.get('satisfied')}", flush=True)

    return True


def main():
    # Must init globals in main process too
    import pyqreach
    pyqreach.initializeTransitionSystem()

    results = []
    for stem, n, file_index, num_qubits in CASES:
        family = "pe" if stem.startswith("pe") else "qft"
        print(f"\n--- {stem} ({family}) ---", flush=True)

        ctx = mp.get_context("spawn")
        # We use a subprocess approach
        proc = ctx.Process(target=worker, args=(stem, n, file_index, num_qubits, family))
        proc.start()
        proc.join(300)
        if proc.is_alive():
            proc.terminate()
            proc.join()
            results.append((stem, "TIMEOUT"))
            print(f"  -> {stem}: TIMEOUT", flush=True)
        else:
            code = proc.exitcode
            if code == 0:
                results.append((stem, "OK"))
                print(f"  -> {stem}: OK (no crash)", flush=True)
            elif code == -11 or code == 139:
                results.append((stem, "SIGSEGV"))
                print(f"  -> {stem}: SIGSEGV (crash!)", flush=True)
            elif code == -10 or code == 138:
                results.append((stem, "SIGBUS"))
                print(f"  -> {stem}: SIGBUS (crash!)", flush=True)
            else:
                results.append((stem, f"EXIT={code}"))
                print(f"  -> {stem}: EXIT={code}", flush=True)

    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    for stem, status in results:
        marker = "<<< CRASH" if status in ("SIGSEGV", "SIGBUS") else ""
        print(f"  {stem:10s}  {status:12s}  {marker}")


if __name__ == "__main__":
    main()
