from collections import defaultdict
import numpy as np
from booldog.utils.misc import ensure_ndarray
from booldog.utils import boolean_normal_forms
import logging

logger = logging.getLogger(__name__)


##############################
#     SQUAD/INTERACTIONS     #
##############################


def _conform_interactions(interactions, activator_value, inhibitor_value):

    if (activator_value != 1) or (inhibitor_value != -1):
        sign_translate = {
            activator_value: 1,
            inhibitor_value: -1,
        }
        for s in interactions:
            for t in interactions[s]:
                symbol = interactions[s][t]
                if symbol in sign_translate:
                    interactions[s][t] = sign_translate[symbol]
                else:
                    logger.warning(
                        f'Issue with edge {s} --> {t}: '
                        f'"{symbol}" not recognized as activator or inhibitor, '
                        f'perhaps you need to define the '
                        f'"inhibitor_value" ({inhibitor_value}) or '
                        f'"activator_value" ({activator_value}) in keyword arguments. '
                    )


class SquadInteractions:

    def __init__(self, interactions, activator_value=1, inhibitor_value=-1, **_):
        ''' '''

        _conform_interactions(interactions, activator_value, inhibitor_value)

        print(interactions)

        self.nodes = tuple(sorted(interactions.keys())) # tuple (i.e. immutable)
        self.n = len(self.nodes)
        self.index = {i:node for node, i in enumerate(self.nodes)}

        funcs = self.squad_update_funcs(interactions)
        self.primes = boolean_normal_forms.functions2primes(functions=funcs)


    def _make_reverse_and_complete(self, interactions):
        reverse_interactions = defaultdict(dict)

        for source, d in interactions.items():
            for target, sign in d.items():
                reverse_interactions[target][source] = sign

        # for node in self.nodes:
        #     regulators = interactions[node]
        #     if len(regulators) == 0:
        #         interactions[node][node] = '+'
        return reverse_interactions

    def interactions_to_matrices(self, reverse_interactions):
        '''
        Only if graph is of type threshold (i.e. SQUAD) does this make sense.
        Create logic matrices
        TODO: these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to
        use e.g. dok_matrix
        '''
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors

        for target, d in reverse_interactions.items():
            for source, sign in d.items():
                if sign == 1:
                    Act[self.index[target], self.index[source]] = 1
                elif sign == -1:
                    Inh[self.index[target], self.index[source]] = 1
                else:
                    print("Warning: Issue with edge: ", source, target)
        return ensure_ndarray(Act), ensure_ndarray(Inh)

    def squad_update_funcs(self, interactions):
        # TODO remove dep on act and inh

        funcs = {}

        # interactions are [target][source]
        reverse_interactions = self._make_reverse_and_complete(interactions)

        Act, Inh = self.interactions_to_matrices(reverse_interactions)

        for node, d in reverse_interactions.items():
            args = list(d.keys())

            def func(*func_input, node=node, args=args):
                if len(func_input) != len(args):
                    print("an issue happened #1")
                    print(node, func.depends, args, func_input)

                state = np.ones(self.n)
                for other_node, other_node_state in zip(args, func_input):
                    state[self.index[other_node]] = other_node_state

                inh = Inh[self.index[node],:].dot(state) > 0
                act = Act[self.index[node],:].dot(state) > 0
                node_state = act * (1-inh)
                return node_state
            func.node = node
            func.depends = args

            funcs[node] = func

        return funcs