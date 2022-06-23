##############################
#     SQUAD/INTERACTIONS     #
##############################

class SquadInteractions:

    def __init__(self, data):
        self.nodes = tuple(sorted(data.keys())) # tuple (i.e. immutable)
        self.n = len(self.nodes)
        self.index = {i:node for node, i in enumerate(self.nodes)}

        funcs = self.squad_update_funcs(data)
        self.primes = boolean_normal_forms.functions2primes(funcs)


    def _make_reverse_and_complete(self, data):
        interactions = defaultdict(dict)

        for source, d in data.items():
            for target, sign in d.items():
                interactions[target][source] = sign

        # for node in self.nodes:
        #     regulators = interactions[node]
        #     if len(regulators) == 0:
        #         interactions[node][node] = '+'
        return interactions

    def interactions_to_matrices(self, data):
        '''
        Only if graph is of type threshold (i.e. SQUAD) does this make sense.
        Create logic matrices
        TODO: these can be made sparse matrices without much effort
        see https://docs.scipy.org/doc/scipy/reference/sparse.html to
        use e.g. dok_matrix
        '''
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors

        for target, d in data.items():
            for source, sign in d.items():
                if sign == "+":
                    Act[self.index[target], self.index[source]] = 1
                elif sign == "-":
                    Inh[self.index[target], self.index[source]] = 1
                else:
                    print("Warning: Issue with edge: ", source, target)
        return ensure_ndarray(Act), ensure_ndarray(Inh)

    def squad_update_funcs(self, interactions):
        # TODO remove dep on act and inh

        funcs = {}

        # interactions are [target][source]
        interactions = self._make_reverse_and_complete(interactions)

        Act, Inh = self.interactions_to_matrices(interactions)

        for node, d in interactions.items():
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