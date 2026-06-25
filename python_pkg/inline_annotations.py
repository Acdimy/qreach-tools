"""Lightweight Qiskit helpers for QReach inline annotations.

The main entry point is QReachCircuit, a small wrapper around Qiskit's
QuantumCircuit.  It delegates normal circuit operations to the wrapped circuit
and records user marks as labelled barrier instructions.  parse_qiskit.py treats
those barriers as semantic no-ops and records the current transition-system
locations for later quantum/state annotations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from qiskit import QuantumCircuit


MARK_LABEL_PREFIX = "qreach:mark:"


def mark(circuit: QuantumCircuit | "QReachCircuit", name: str) -> QuantumCircuit | "QReachCircuit":
    """Insert a QReach mark at the current position in a circuit.

    The mark is represented as a labelled barrier over all qubits.  It is a
    no-op for QReach's transition-system semantics, but parse_qiskit_cir records
    the current locations under `name`.
    """
    if isinstance(circuit, QReachCircuit):
        circuit.mark(name)
        return circuit
    circuit.barrier(*list(circuit.qubits), label=f"{MARK_LABEL_PREFIX}{name}")
    return circuit


@dataclass
class QReachCircuit:
    """A thin wrapper that supports normal Qiskit syntax plus .mark(name).

    Unknown attributes and methods are delegated to the underlying
    QuantumCircuit, so most existing Qiskit code can be used unchanged after
    constructing QReachCircuit instead of QuantumCircuit.  Use `.circuit` or
    `.unwrap()` when an API requires a real QuantumCircuit instance.
    """

    circuit: QuantumCircuit
    declared_marks: list[str] = field(default_factory=list)

    def __init__(self, *args: Any, **kwargs: Any):
        if len(args) == 1 and isinstance(args[0], QuantumCircuit) and not kwargs:
            self.circuit = args[0]
        else:
            self.circuit = QuantumCircuit(*args, **kwargs)
        self.declared_marks = []

    def __getattr__(self, name: str):
        return getattr(self.circuit, name)

    def __len__(self) -> int:
        return len(self.circuit)

    def __iter__(self):
        return iter(self.circuit)

    def __getitem__(self, item):
        return self.circuit[item]

    def __repr__(self) -> str:
        return f"QReachCircuit({self.circuit!r})"

    @classmethod
    def from_qasm_file(cls, path: str) -> "QReachCircuit":
        """Load a normal OpenQASM 2 file and wrap it as a QReachCircuit.

        This mirrors ``QuantumCircuit.from_qasm_file``.  The QASM file is assumed
        to be ordinary OpenQASM without QReach-specific inline annotations; users
        may still add ``.mark(name)`` calls to the returned wrapper before parsing.
        """
        return cls(QuantumCircuit.from_qasm_file(path))

    @classmethod
    def from_qasm_str(cls, qasm_str: str) -> "QReachCircuit":
        """Load an OpenQASM 2 string and wrap it as a QReachCircuit."""
        return cls(QuantumCircuit.from_qasm_str(qasm_str))

    def mark(self, name: str) -> "QReachCircuit":
        mark(self.circuit, name)
        self.declared_marks.append(name)
        return self

    def unwrap(self) -> QuantumCircuit:
        return self.circuit


def is_mark_operation(operation) -> bool:
    label = getattr(operation, "label", None)
    return operation.name == "barrier" and isinstance(label, str) and label.startswith(MARK_LABEL_PREFIX)


def mark_name(operation) -> str:
    label = getattr(operation, "label", "")
    if not isinstance(label, str) or not label.startswith(MARK_LABEL_PREFIX):
        raise ValueError("Operation is not a QReach mark")
    return label[len(MARK_LABEL_PREFIX):]
