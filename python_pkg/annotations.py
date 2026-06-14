"""Built-in annotation predicates for QReach transition systems.

This module provides a small registry-based API for turning common transition
system predicates into NuSMV/qCTL atomic proposition labels.  It intentionally
uses the public pyqreach/qctl-style interfaces so it can work with both the
explicit TransitionSystem and, where supported, SymTS.
"""

from __future__ import annotations

from dataclasses import dataclass
from fnmatch import fnmatch
from typing import Callable, Iterable, Mapping


Predicate = Callable[[object, int], bool]


@dataclass(frozen=True)
class AnnotationKeyword:
    """A named transition-system predicate that can label locations."""

    name: str
    predicate: Predicate
    description: str = ""


def _location_ids(ts) -> list[int]:
    if hasattr(ts, "getLocationIDs"):
        return list(ts.getLocationIDs())
    return list(range(ts.getLocationNum()))


def _get_identifier(ts, loc: int) -> str:
    if hasattr(ts, "getIdentifier"):
        return ts.getIdentifier(loc)
    return ts.Locations[loc].getIdentifier()


def _get_post_locations(ts, loc: int) -> list[int]:
    if hasattr(ts, "getPostLocations"):
        return list(ts.getPostLocations(loc))
    return list(ts.Locations[loc].postLocations)


def _relation_name(ts, src: int, dst: int) -> str:
    try:
        return str(ts.getRelationName(src, dst))
    except Exception:
        return ""


def _satisfy_bit(ts, loc: int, idxs: list[int], vals: list[int]) -> list[str]:
    if hasattr(ts, "satisfyBit"):
        return list(ts.satisfyBit(loc, idxs, vals))
    return list(ts.Locations[loc].satisfyBit(idxs, vals))


def _is_reached(ts, loc: int) -> bool:
    """Return whether fixed-point computation found a non-zero state at loc.

    This follows qctl.tsLabellingDefault: a location is reached when
    ts.printDims(loc)[1] > 0.  For the explicit qts_naive::TransitionSystem
    this is the canonical existing workflow-level predicate.
    """
    try:
        return ts.printDims(loc)[1] > 0
    except Exception:
        # TODO: SymTS has locationHasNonZeroAnnotation(loc), but the public
        # workflow should converge on one transition-system-independent reached
        # predicate instead of branching here.
        if hasattr(ts, "locationHasNonZeroAnnotation"):
            try:
                return bool(ts.locationHasNonZeroAnnotation(loc))
            except Exception:
                pass
        return False


def _is_leaf(ts, loc: int) -> bool:
    try:
        return bool(ts.isLeafLoc(loc))
    except Exception:
        return len(_get_post_locations(ts, loc)) == 0


def _is_init(ts, loc: int) -> bool:
    return loc == ts.getInitLocation()


def _is_loop(ts, loc: int) -> bool:
    # Current parser identifiers use "W" for while-loop-generated locations.
    # TODO: replace this heuristic with structured parser metadata once
    # parse_qiskit_cir returns ParseResult/marker information.
    return "W" in _get_identifier(ts, loc)


def _is_deadend(ts, loc: int) -> bool:
    # This is intentionally graph-topological.  It is not the same as
    # unreachable/non-reached.
    return len(_get_post_locations(ts, loc)) == 0


def _is_measured(ts, loc: int) -> bool:
    # A measured location is currently inferred from an incoming meas0/meas1
    # transition.  TODO: have parse_qiskit.py attach a structured measurement
    # marker instead of recovering it from relation names.
    for src in _location_ids(ts):
        if loc in _get_post_locations(ts, src):
            name = _relation_name(ts, src, loc).lower()
            if "meas0" in name or "meas1" in name:
                return True
    return False


def _is_branch(ts, loc: int) -> bool:
    # Coarse graph predicate: any location with multiple outgoing targets.
    # TODO: distinguish measurement branches from if/while/switch control-flow
    # branches when parser metadata is available.
    return len(set(_get_post_locations(ts, loc))) > 1


class AnnotationRegistry:
    """Registry for transition-system annotation keywords.

    Registered predicates are applied as labels on locations.  The labels then
    become atomic propositions in qctl.dict2SMV/modelChecking.
    """

    def __init__(self):
        self._keywords: dict[str, AnnotationKeyword] = {}

    @classmethod
    def builtins(cls) -> "AnnotationRegistry":
        registry = cls()
        registry.register("reached", _is_reached, "location has non-zero fixed-point annotation")
        registry.register("valid", _is_reached, "alias for reached, kept for existing scripts")
        registry.register("leaf", _is_leaf, "terminal/leaf transition-system location")
        registry.register("init", _is_init, "initial transition-system location")
        registry.register("loop", _is_loop, "location generated inside a while-loop region")
        registry.register("deadend", _is_deadend, "location with no outgoing graph edges")
        registry.register("measured", _is_measured, "location reached by a measurement transition")
        registry.register("branch", _is_branch, "location with multiple outgoing graph edges")
        return registry

    def register(
        self,
        name: str,
        predicate: Predicate,
        description: str = "",
        *,
        overwrite: bool = False,
    ) -> None:
        if not overwrite and name in self._keywords:
            raise ValueError(f"Annotation keyword already registered: {name}")
        self._keywords[name] = AnnotationKeyword(name, predicate, description)

    def register_identifier(self, name: str, pattern: str, *, overwrite: bool = False) -> None:
        """Register a keyword that labels locations whose identifier matches a glob."""
        self.register(
            name,
            lambda ts, loc: fnmatch(_get_identifier(ts, loc), pattern),
            f"identifier matches {pattern!r}",
            overwrite=overwrite,
        )

    def register_classical(
        self,
        name: str,
        bit_patterns: str | Iterable[str],
        *,
        bit_indices: Iterable[int] | None = None,
        overwrite: bool = False,
    ) -> None:
        """Register a keyword for classical-bit valuations.

        bit_patterns follows the existing tsLabellingClRegList convention: each
        string is a bit pattern such as "10", and the keyword is true if any
        pattern is satisfied.  By default the patterns are matched against the
        first len(pattern) classical bits.
        """
        if isinstance(bit_patterns, str):
            patterns = [bit_patterns]
        else:
            patterns = list(bit_patterns)

        def predicate(ts, loc: int) -> bool:
            for pattern in patterns:
                vals = [int(bit) for bit in pattern]
                idxs = list(bit_indices) if bit_indices is not None else list(range(len(vals)))
                if len(_satisfy_bit(ts, loc, idxs, vals)) != 0:
                    return True
            return False

        self.register(name, predicate, f"classical bit pattern in {patterns!r}", overwrite=overwrite)

    def names(self) -> list[str]:
        return sorted(self._keywords.keys())

    def get(self, name: str) -> AnnotationKeyword:
        try:
            return self._keywords[name]
        except KeyError as exc:
            raise KeyError(f"Unknown annotation keyword: {name}") from exc

    def apply(
        self,
        ts,
        names: Iterable[str] | Mapping[str, Predicate] | None = None,
        *,
        loc_list: Iterable[int] | None = None,
    ) -> dict[str, list[int]]:
        """Apply annotation keywords to a transition system.

        Returns a mapping from label name to the list of locations that were
        labelled.  `names=None` applies all registered keywords.  A mapping may
        be passed to apply ad-hoc predicates without registering them first.
        """
        locs = list(loc_list) if loc_list is not None else _location_ids(ts)
        predicates: dict[str, Predicate]
        if names is None:
            predicates = {name: keyword.predicate for name, keyword in self._keywords.items()}
        elif isinstance(names, Mapping):
            predicates = dict(names)
        else:
            predicates = {name: self.get(name).predicate for name in names}

        labelled: dict[str, list[int]] = {name: [] for name in predicates}
        for name, predicate in predicates.items():
            for loc in locs:
                if predicate(ts, loc):
                    ts.setLabel(loc, name)
                    labelled[name].append(loc)
        return labelled


_DEFAULT_REGISTRY = AnnotationRegistry.builtins()


def default_registry() -> AnnotationRegistry:
    """Return the process-wide default annotation registry."""
    return _DEFAULT_REGISTRY


def annotate(ts, names: Iterable[str] | Mapping[str, Predicate] | None = None, *, loc_list=None):
    """Apply built-in or ad-hoc annotations using the default registry."""
    return _DEFAULT_REGISTRY.apply(ts, names, loc_list=loc_list)


def annotate_where(ts, name: str, predicate: Predicate, *, loc_list=None) -> dict[str, list[int]]:
    """Apply one ad-hoc annotation predicate without registering it globally."""
    return _DEFAULT_REGISTRY.apply(ts, {name: predicate}, loc_list=loc_list)


def annotate_identifier(ts, name: str, pattern: str, *, loc_list=None) -> dict[str, list[int]]:
    """Apply an identifier-glob annotation without registering it globally."""
    return annotate_where(ts, name, lambda target_ts, loc: fnmatch(_get_identifier(target_ts, loc), pattern), loc_list=loc_list)


def annotate_classical(
    ts,
    name: str,
    bit_patterns: str | Iterable[str],
    *,
    bit_indices: Iterable[int] | None = None,
    loc_list=None,
) -> dict[str, list[int]]:
    """Apply a classical-bit-pattern annotation without registering it globally."""
    if isinstance(bit_patterns, str):
        patterns = [bit_patterns]
    else:
        patterns = list(bit_patterns)

    def predicate(target_ts, loc: int) -> bool:
        for pattern in patterns:
            vals = [int(bit) for bit in pattern]
            idxs = list(bit_indices) if bit_indices is not None else list(range(len(vals)))
            if len(_satisfy_bit(target_ts, loc, idxs, vals)) != 0:
                return True
        return False

    return annotate_where(ts, name, predicate, loc_list=loc_list)
