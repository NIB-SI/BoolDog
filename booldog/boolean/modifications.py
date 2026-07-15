'''Module containing functions to manipulate Boolean networks.
It is not intended to be used directly, but rather as a mixin
for other classes.
'''
import logging

from pyboolnet.prime_implicants import find_successors, create_variables, remove_variables

from booldog.utils import ExtendedEnum
from booldog.classes import BoolDogNode

logger = logging.getLogger(__name__)

class ModificationTypes(ExtendedEnum):
    '''Types of modifications to the Boolean network.

    :meta private:
    '''
    ADD = "add_node"
    REMOVE = "remove_node"
    UPDATE = "update"


class Modification():
    '''Represents a single modification made (or to be made) to a Boolean
    network. Used both as the audit-trail record appended to
    `BoolDogModel.modifications` by `_track_network_modifications`, and as
    the input format accepted by `modify_network`.
    '''

    def __init__(self, modification_type, node_id, rule=None):
        self.type = modification_type
        '''ModificationTypes : the type of modification (add, remove, or update).'''
        self.node_id = node_id
        '''str or list of str : identifier(s) of the node(s) affected.'''
        self.rule = rule
        '''str or None : the rule (bnet form) associated with the modification, if any.'''

    def __repr__(self):
        return f"Modification(type={self.type}, node_id={self.node_id}, rule={self.rule})"


class BooleanNetworkModificationMixin():
    '''Mixin class to modify a Boolean network.
    This class is not intended to be used directly, but rather as a mixin.
    '''
    def _track_network_modifications(self, modification_type, node_id, rule=None):
        '''
        Append a `Modification` record to `self.modifications`, the model's
        audit trail of changes made via `add_node`/`remove_node(s)`/
        `update_node`.

        Parameters
        ----------
        modification_type : ModificationTypes
            Type of modification performed.
        node_id : str or list of str
            Identifier(s) of the node(s) affected.
        rule : str, optional
            The rule (bnet form) associated with the modification, if any.
        '''
        self.modifications.append(
            Modification(modification_type, node_id, rule=rule))

    # function to update the object attributes
    def _update_model_object(self, uncache_primes=True):
        '''Refresh the model's derived bookkeeping after a node addition,
        removal, or rule update.

        Recomputes `self.node_ids`/`self.index` (via
        `BoolDogModel._set_node_ids_and_index`) and, unless the caller has
        already updated `self.primes` in place (e.g. via
        `pyboolnet.prime_implicants.create_variables`/`remove_variables`),
        invalidates the cached primes so they are recomputed from the
        node rules the next time `primes` is accessed.

        Parameters
        ----------
        uncache_primes : bool, optional
            Whether to invalidate the primes cache (default True). Pass
            False when the primes have already been updated in place by
            the caller.
        '''

        self._set_node_ids_and_index() # Mixin function

        if uncache_primes:
            # next time primes is accessed, it will be recalculated from the rules
            self._primes_cached = False

    def modify_network(self, modifications):
        '''
        Modify the network according to the given modifications.

        Parameters
        ----------
        modifications : list of Modification
            List of `Modification` objects, applied in the order given.

        Returns
        ----------
        None
            The network is modified in place.

        Notes
        -----

        Each `Modification` object has the following attributes:

        - type : ModificationTypes
          Type of modification. One of `ModificationTypes.ADD`
          ("add_node"), `ModificationTypes.REMOVE` ("remove_node"), or
          `ModificationTypes.UPDATE` ("update").
        - node_id : str or list
          Identifier(s) of node(s) to modify.
        - rule : str or None
          Rule to define the update function of `node_id`, in bnet form.
          All nodes in `rule` need to be defined in Network.

        The modifications are applied in the order they are given, by
        dispatching each one to `add_node`, `remove_nodes`, or (for
        "update") the equivalent of `update_node`.
        A node cannot be removed if other nodes depend on it (i.e. it occurs in
        their update logic). To remove such a node, either also remove all of
        its dependants, or first update the logic rule of its dependants to
        remove dependency.

        This is a wrapper for pyboolnet.prime_implicants.create_variables and
        pyboolnet.prime_implicants.remove_variables.
        '''

        if not isinstance(modifications, list):
            raise ValueError(
                "Modifications must be a list of Modification objects.")
        if not all(isinstance(x, Modification) for x in modifications):
            raise ValueError("All modifications must be Modification objects.")

        for modification in modifications:
            match modification.type:
                case ModificationTypes.ADD:
                    self.add_node(modification.node_id, modification.rule)
                case ModificationTypes.REMOVE:
                    self.remove_nodes(modification.node_id)
                case ModificationTypes.UPDATE:
                    self.update_node(modification.node_id, modification.rule)
                case _:
                    raise ValueError(
                        f"Modification type has to be one of {', '.join(ModificationTypes.values())}, not {modification.type}."
                    )

    def remove_node(self, node_id):
        '''
        Removes `node_id` from the network. Convenience wrapper around
        `remove_nodes` for a single node.

        Parameters
        ----------
        node_id : str
            Identifier of the node to remove.

        Returns
        ----------
        None

        Notes
        -----
        A node cannot be removed if other nodes depend on it (i.e. it occurs in
        their update logic). To remove such a node, either also remove all of
        its dependants, or first update the logic rule of its dependants to
        remove dependency.

        This is a wrapper for pyboolnet.prime_implicants.remove_variables.
        '''

        self.remove_nodes(node_id)

    def _test_node_removabilty(self, node_ids):
        '''Check that removing all of `node_ids` from the network at once is
        valid, i.e. that no node outside `node_ids` depends on (is a
        successor of) any node in `node_ids`.

        Parameters
        ----------
        node_ids : list of str
            Identifiers of the nodes to be removed together.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If any node in `node_ids` has a dependant (successor, via
            `pyboolnet.prime_implicants.find_successors`) that is not
            itself also being removed. The message lists, for each such
            node, which of its dependants block the removal.
        '''

        hit = {
            node_id: [
                x for x in find_successors(primes=self.primes, sources=[node_id])
                if x not in node_ids
            ]
            for node_id in node_ids
        }
        s = ''
        for node_id in hit:
            if hit[node_id]:  # if there are dependants
                s += f'Cannot remove a node that has dependents. To remove '\
                f'"{node_id}", you need to remove its dependants as well, or '\
                f'remove their dependency on "{node_id}" by updating their rules. '\
                f'Dependent(s) for "{node_id}" are: {", ".join(hit[node_id])}.\n'
                # TODO example

        if s:
            raise ValueError(s)

    def remove_nodes(self, node_ids):
        '''
        Removes all nodes in `node_ids` from the network.

        Parameters
        ----------
        node_ids : str or list of str
            Identifier, or list of identifiers, of nodes to remove. A
            single string is treated as one node identifier.

        Returns
        ----------
        None

        Raises
        ------
        ValueError
            If any of `node_ids` is not present in the network, or if
            removing all of `node_ids` together would leave a dependant
            node whose rule still refers to a removed node (see
            `_test_node_removabilty`).

        Notes
        -----
        A node cannot be removed if other nodes depend on it (i.e. it occurs in
        their update logic). To remove such a node, either also remove all of
        its dependants, or first update the logic rule of its dependants to
        remove dependency.

        This is a wrapper for pyboolnet.prime_implicants.remove_variables.
        '''
        if isinstance(node_ids, str):
            node_ids = [node_ids]

        if not all(node_id in self.node_ids for node_id in node_ids):
            raise ValueError(
                f"All nodes to remove must be present in the network. Nodes not found: {', '.join([node_id for node_id in node_ids if node_id not in self.node_ids])}."
            )

        self._test_node_removabilty(node_ids)

        remove_variables(self.primes, node_ids)

        # remove the node from the model
        for node_id in node_ids:
            del self.nodes[node_id]
        self._update_model_object(uncache_primes=False) # remove_variables already updated primes

        self._track_network_modifications(ModificationTypes.REMOVE, node_ids)

    def add_node(self, node_id, rule, name=None):
        '''
        Add a new node to Network.

        Parameters
        ----------
        node_id : str
            Identifier of the new node. Must not already be present in
            the network.
        rule : str
            Rule to define update of `node_id`, in bnet form. All nodes in
            `rule` need to be defined in Network.
        name : str, optional
            Not currently used by this method (the created `BoolDogNode`
            is constructed without a `name`, so it falls back to
            `node_id`; see `BoolDogNode.__post_init__`).

        Returns
        ----------
        None

        Raises
        ------
        ValueError
            If `node_id` is already present in the network.

        Notes
        -----
        This is a wrapper for pyboolnet.prime_implicants.create_variables.
        '''
        if node_id in self.nodes:
            raise ValueError(f"{node_id} is already present in Network.")

        self.nodes[node_id] = BoolDogNode(identifier=node_id, rule=rule)

        create_variables(self.primes, {node_id: rule})

        self._update_model_object(uncache_primes=False) # create_variables already updated primes

        self._track_network_modifications(ModificationTypes.ADD, node_id, rule=rule)



    def update_node(self,
                    node_id,
                    rule,
                    modification_type=ModificationTypes.UPDATE):
        '''
        Update (overwrite) the logic rule defining the update of `node_id`.
        `node_id` must already exist in the network; use `add_node` to add
        a new node instead.

        Parameters
        ----------
        node_id : str
            Identifier of the (existing) node whose rule is to be updated.
        rule : str
            New rule to define update of `node_id`, in bnet form. All nodes
            in `rule` need to be defined in Network.
        modification_type : ModificationTypes, optional
            The modification type recorded in the audit trail
            (`self.modifications`) for this change. Defaults to
            `ModificationTypes.UPDATE`.

        Returns
        ----------
        None

        Raises
        ------
        ValueError
            If `node_id` is not already present in the network.

        Notes
        -----
        If `rule` is identical to the node's current rule, a warning is
        logged and no update (and no audit-trail entry) is made.

        This is a wrapper for pyboolnet.prime_implicants.create_variables.
        '''
        if node_id not in self.nodes:
            raise ValueError(f"{node_id} is not present in Network. To add a new node, use add_node().")

        if self.nodes[node_id].rule == rule:
            logger.warning("Rule for node '%s' is already '%s'. No update performed.", node_id, rule)
            return

        create_variables(self.primes, {node_id: rule})

        self.nodes[node_id].rule = rule

        self._update_model_object(uncache_primes=False) # create_variables already updated primes

        self._track_network_modifications(modification_type, node_id, rule=rule)
