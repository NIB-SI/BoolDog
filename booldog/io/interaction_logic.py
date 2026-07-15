''' Convert interactions into Boolean functions
'''

import logging
from collections import defaultdict
from typing import Dict, Mapping

logger = logging.getLogger(__name__)

# ---- Public API -------------------------------------------------------------

def interactions2rules(
    interactions: tuple[int | str, int | str, int | str],
    *,
    logic: "LogicBuilder" = None,
    activator_symbol: int | str = 1,
    inhibitor_symbol: int | str = -1,
) -> Dict[str, str]:
    '''Convert a list of signed interactions into Boolean update rules.

    Interactions are first grouped by target node (see
    `_normalise_and_collect_regulators`), then each target's Boolean rule
    string (bnet-format, e.g. "A | B & !C") is built from its regulators by
    `logic`.

    Parameters
    ----------
    interactions : list
        A list of interactions in the network. Each interaction should be a tuple of
        (source, target, sign), where source and target are node identifiers and
        sign is either activator_symbol or inhibitor_symbol.

    activator_symbol : int, optional
        The value representing activation in the network. Default is 1.

    inhibitor_symbol : int, optional
        The value representing inhibition in the network. Default is -1.

    logic : LogicBuilder, optional
        An optional logic builder to use for constructing the update functions.
        If not provided, the default is `SquadLogic`, which implements the SQUAD logic:
        A node is active iff::

            (any activator is active) AND (no inhibitor is active)


    Returns
    -------
    dict
        A dictionary mapping each target node identifier to its Boolean
        update rule, as a bnet-format rule string (not a prime-implicant
        object) produced by ``logic.build``.

    '''

    logic = logic or SquadLogic()

    regulators_per_target = _normalise_and_collect_regulators(
        interactions,
        activator_symbol=activator_symbol,
        inhibitor_symbol=inhibitor_symbol,
    )

    rules = {}
    for node, regulators in regulators_per_target.items():
        rules[node] = logic.build(node, regulators)

    return rules

# ---- Interaction preprocessing ----------------------------------------------

def _normalise_and_collect_regulators(interactions, activator_symbol, inhibitor_symbol):
    """Normalize signs and orient interactions as target -> regulators.

    Parameters
    ----------
    interactions : list of tuple
        Interactions as (source, target, sign) tuples.
    activator_symbol : int or str
        The value of `sign` that represents activation; normalized to ``1``.
    inhibitor_symbol : int or str
        The value of `sign` that represents inhibition; normalized to ``-1``.

    Returns
    -------
    dict
        A dictionary mapping each target node identifier to a dict of its
        regulators, i.e. ``{target: {regulator: 1 or -1, ...}, ...}``.
        Interactions whose `sign` matches neither `activator_symbol` nor
        `inhibitor_symbol` are dropped (a warning is logged) rather than
        included in the result. If the same (source, target) pair appears
        more than once, only the last occurrence is kept.
    """
    translate = {
        activator_symbol: 1,
        inhibitor_symbol: -1,
    }

    regulators_per_target = defaultdict(dict)

    for source, target, sign in interactions:
        norm = translate.get(sign)
        if norm is None:
            logger.warning(
                "Ignoring edge %s → %s with unrecognized sign %r",
                source,
                target,
                sign,
            )
            continue
        regulators_per_target[target][source] = norm

    return dict(regulators_per_target)

# ---- Logic interface --------------------------------------------------------

class LogicBuilder:
    """Abstract base class for Boolean logic builders.

    Subclasses implement `build` to turn a node's set of signed regulators
    into a Boolean rule string (bnet format), pluggable into
    `interactions2rules`.
    """

    def build(self, node: str, regulators: Mapping[str, int]) -> str:
        """Return an update rule for `node`.

        Parameters
        ----------
        node : str
            Identifier of the target node the rule is being built for.
        regulators : Mapping[str, int]
            Mapping of regulator node identifier to normalized sign
            (``1`` for activator, ``-1`` for inhibitor).

        Returns
        -------
        str
            A Boolean rule string (bnet format) for `node`.
        """
        raise NotImplementedError

class SquadLogic(LogicBuilder):
    """SQUAD logic:

    A node is active iff::

        (any activator is active) AND (no inhibitor is active)
    """

    def build(self, node: str, regulators: Mapping[str, int]) -> str:
        """Build the SQUAD-logic rule for `node` from its regulators.

        Parameters
        ----------
        node : str
            Identifier of the target node the rule is being built for
            (unused other than for interface compatibility with
            `LogicBuilder.build`).
        regulators : Mapping[str, int]
            Mapping of regulator node identifier to normalized sign
            (``1`` for activator, ``-1`` for inhibitor).

        Returns
        -------
        str
            A Boolean rule string (bnet format): ``"0"`` if there are no
            regulators; the activators OR-ed together if there are no
            inhibitors; the negated inhibitors AND-ed together if there are
            no activators; otherwise ``"(activators OR-ed) & (negated
            inhibitors AND-ed)"``.
        """
        activators = [r for r, s in regulators.items() if s == 1]
        inhibitors = [r for r, s in regulators.items() if s == -1]

        if not activators and not inhibitors:
            return "0"
        elif not inhibitors:
            return " | ".join(activators)
        elif not activators:
            return " & ".join(f"!{r}" for r in inhibitors)
        else:
            return f"({' | '.join(activators)}) & ({' & '.join(f'!{r}' for r in inhibitors)})"
