import pyqreach
from qiskit import QuantumCircuit
from parse_qiskit import parse_qiskit_cir_lazy


def _dims(ts, loc):
    return ts.printDims(loc)


def test_deterministic_measurement_prunes_zero_outcomes():
    pyqreach.initializeTransitionSystem()
    qc = QuantumCircuit(2, 2)
    qc.measure(0, 0)
    qc.measure(1, 1)

    ts = pyqreach.TransitionSystem()
    result = parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state="00", return_metadata=True)

    assert result.lazy is True
    assert len(result.result_locations) == 1
    assert len(result.lazy_pruned_locations) == 2
    assert ts.getLocationNum() == 5

    for loc in result.lazy_pruned_locations:
        assert _dims(ts, loc)[1] == 0
        assert ts.isLeafLoc(loc)

    reached_leaf = result.result_locations[0]
    assert _dims(ts, reached_leaf)[1] > 0
    assert ts.isLeafLoc(reached_leaf)


def test_hadamard_measurement_keeps_both_outcomes():
    pyqreach.initializeTransitionSystem()
    qc = QuantumCircuit(2, 1)
    qc.h(0)
    qc.measure(0, 0)

    ts = pyqreach.TransitionSystem()
    result = parse_qiskit_cir_lazy(qc, qc.num_qubits, ts, initial_state="00", return_metadata=True)

    assert result.lazy is True
    assert result.lazy_pruned_locations == []
    assert len(result.result_locations) == 2
    for loc in result.result_locations:
        assert _dims(ts, loc)[1] > 0


if __name__ == "__main__":
    test_deterministic_measurement_prunes_zero_outcomes()
    test_hadamard_measurement_keeps_both_outcomes()
    print("All lazy measurement tests PASSED")
