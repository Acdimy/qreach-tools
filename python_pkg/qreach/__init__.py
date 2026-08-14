"""QReach Python workflow package.

Provides the Qiskit -> TransitionSystem lowering (``parse_qiskit``) and the
model-checking workflow utilities (``qctl``), plus annotation helpers and
quantum error-correcting circuit builders. Import submodules directly, e.g.::

    from qreach.qctl import quantum_state, span_qops, modelChecking
    from qreach.parse_qiskit import parse_qiskit_cir, parse_qiskit_cir_lazy
"""

__all__ = ["symbolic_available"]


def symbolic_available() -> bool:
    """Return whether the symbolic SymTS backend is available.

    The symbolic transition-system path is CFLOBDD-only; under the LimTDD
    backend the SymTS bindings are compiled out, so this returns False and the
    symbolic regression tests should skip.
    """
    try:
        import pyqreach
    except ImportError:
        return False
    return hasattr(pyqreach, "SymTS")
