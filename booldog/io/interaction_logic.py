''' Convert interactions into Boolean functions
'''

import logging
from collections import defaultdict
from typing import Callable, Dict, Mapping

logger = logging.getLogger(__name__)

# ---- Public API -------------------------------------------------------------

def interactions2rules(
    interactions: tuple[int | str, int | str, int | str],
    *,
    logic: "LogicBuilder" = None,
    activator_symbol: int | str = 1,
    inhibitor_symbol: int | str = -1,
) -> Dict[str, Callable]:
    '''Convert interactions directly to prime implicants.

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
        A node is active iff:
            (any activator is active) AND (no inhibitor is active)


    Returns
    -------
    dict
        A dictionary mapping node identifiers to their corresponding prime implicants.

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
    """Normalize signs and orient interactions as target → regulators."""
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
    """Abstract base class for Boolean logic builders."""

    def build(self, node: str, regulators: Mapping[str, int]) -> str:
        """Return an update rule for `node`."""
        raise NotImplementedError

class SquadLogic(LogicBuilder):
    """SQUAD logic:

    A node is active iff:
        (any activator is active) AND (no inhibitor is active)
    """

    def build(self, node: str, regulators: Mapping[str, int]) -> str:
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
