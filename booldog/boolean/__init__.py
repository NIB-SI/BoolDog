'''Boolean-network functionality for BoolDogModel.

Exposes :class:`~booldog.boolean.boolean.BooleanNetworkMixin` (read-only
Boolean-network operations: state-space generation, Boolean simulation,
steady states, parent/child lookups) and
:class:`~booldog.boolean.modifications.BooleanNetworkModificationMixin`
(structural edits: `add_node`, `remove_node(s)`, `update_node`,
`modify_network`, and the `Modification`/`ModificationTypes` audit trail).
Both are mixed into :class:`booldog.network.BoolDogModel`.
'''

from .boolean import BooleanNetworkMixin
from .modifications import BooleanNetworkModificationMixin
