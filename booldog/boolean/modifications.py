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
    '''Class to represent a modification to the Boolean network.'''

    def __init__(self, modification_type, node_id, rule=None):
        self.type = modification_type
        self.node_id = node_id
        self.rule = rule

    def __repr__(self):
        return f"Modification(type={self.type}, node_id={self.node_id}, rule={self.rule})"


class BooleanNetworkModificationMixin():
    '''Mixin class to modify a Boolean network.
    This class is not intended to be used directly, but rather as a mixin.
    '''
    def _track_network_modifications(self, modification_type, node_id, rule=None):
        '''
        Keep track of network modifications
        '''
        self.modifications.append(
            Modification(modification_type, node_id, rule=rule))

    # function to update the object attributes
    def _update_model_object(self, uncache_primes=True):
        '''Update the model attributes after a modification.'''

        self._set_node_ids_and_index() # Mixin function

        if uncache_primes:
            # next time primes is accessed, it will be recalculated from the rules
            self._primes_cached = False

    def modify_network(self, modifications):
        '''
        Modify the network according to the given modifications.

        Parameters
        ----------
        modifications : list
            List of Modification objects.

        Returns
        ----------
        None
            The network is modified in place.

        Notes
        -----

        Each Modification object must have the following attributes:
        - modification_type : str
            Type of modification. One of "add_node", "remove_node", "update".
        - node : str or list
            Name of node(s) to modify.
        - rule : str or None
            Rule to define the update function of `node`, in bnet form.
            All nodes in `rule` need to be defined in Network.

        The modifications are applied in the order they are given.
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
            match modification.modification_type:
                case ModificationTypes.ADD:
                    self.add_node(modification.node_id, modification.rule)
                case ModificationTypes.REMOVE:
                    self.remove_nodes(modification.node_id)
                case ModificationTypes.UPDATE:
                    self.update_rule(modification.node_id, modification.rule)
                case _:
                    raise ValueError(
                        f"Modification type has to be one of {', '.join(ModificationTypes.values())}, not {modification.modification_type}."
                    )

    def remove_node(self, node_id):
        '''
        Removes `node` from the network.

        Parameters
        ----------
        node : str
            Names of node to remove

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
        '''Test if a node (set) can be removed from the network.
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
        node_ids : list
            List of identifiers of nodes to remove

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
        node : str
            Node name
        rule : str
            Rule to define update of `node`, in bnet form. All nodes in
            `rule` need to be defined in Network.

        Returns
        ----------
        None

        This is a wrapper for pyboolnet.prime_implicants.create_variables.
        '''
        if node_id in self.nodes:
            raise ValueError("f{node} is already present in Network.")

        self.nodes[node_id] = BoolDogNode(identifier=node_id, rule=rule)

        create_variables(self.primes, {node_id: rule})

        self._update_model_object(uncache_primes=False) # create_variables already updated primes

        self._track_network_modifications(ModificationTypes.ADD, node_id, rule=rule)



    def update_node(self,
                    node_id,
                    rule,
                    modification_type=ModificationTypes.UPDATE):
        '''
        Update (modify) or add the logic rule defining the update of `node`.
        If `node` does not yet exist, it will be added, if it does exist, its
        update logic will be overwritten.

        Parameters
        ----------
        node_id : str
            Node identifier
        rule : str
            New rule to define update of `node`, in bnet form. All nodes in
            `rule` need to be defined in Network.

        Returns
        ----------
        None

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
