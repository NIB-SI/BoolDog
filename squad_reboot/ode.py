
import numpy as np

from .utils import *


from .bool_graph import BooleanGraph




class ODE:

    def __init__(self, graph):

        if not isinstance(graph, BooleanGraph):
            raise TypeError("'graph' argument must be a BooleanGraph object")

        print("Initialising ODE system ... ", end="")
        self.n = len(graph)
        self.boolean_graph = graph

    def event_function(self, t, x, event_t, *args):
        return t - event_t
    event_function.terminal = True






# https://github.com/krumsieklab/Odefy/blob/11d048d550a8f64250ba01f76f5a83048c8be6cf/Odefy-1.20/code/models/CreateCubeCalls.m
class BoolCube(ODE):
    '''
    Transforming Boolean models to continuous models: methodology and application to T-cell receptor signaling

    '''
    def __init__(self, graph, transform, tau=1, n=1, k=1):
        super().__init__(graph)
        transform = transform.lower()
        if transform in ('boolecube', 'boolcube'):
            self.dxdt =  self._get_system(tau, lambda x, n, k: x, n=n, k=k)

        elif transform in ('hill', 'hillcube'):
            self.dxdt =  self._get_system(tau, self.hill, n=n, k=k)

        elif transform in ('normalisedhill', 'normalisedhillcube'):
            self.dxdt =  self._get_system(tau, self.normalised_hill, n=n, k=k)
        print("done. ")

    def _get_system(self, tau, transform_function, off_nodes=[],  **kwargs):

        def dxdt(t, x_array, *args):
            b = self.homologue_b1(transform_function(x_array, **kwargs))
            d = 1/tau * ( b  - x_array)
            for i in off_nodes:
                d[i] = 0
            return d

        return dxdt

    def homologue_b1(self, x_array):

        def F(x_array, x_bool):
            #return {0:1-x_array, 1:x_array}[x_bool]
            return x_array * x_bool + (1 - x_bool)*(1 - x_array)

        B_1 = np.zeros(self.n)
        for node in self.boolean_graph.nodes: # iterate over all nodes
            for prime_dict in self.boolean_graph.primes[node][1]: # sum where B() = 1 part 1
                # in this partly fixed state, B() = 1
                # 1-prime implicant holds fixed node states, we need to iterate
                # over all free node states as well
                 for x_bool in self.boolean_graph.generate_states(fixed=prime_dict): # sum where B() = 1 part 2
                     B_1[self.boolean_graph.index[node]]  += np.product(F(x_array, x_bool))
        return B_1


    def hill(self, x_array, n, k):
        return x_array**n / (x_array**n + k**n)

    def normalised_hill(self, x_array, n, k):
        return self.hill(x_array, n, k) / self.hill(1, n, k)





class ShaoODE(ODE):
    '''
    From Boolean Network Model to Continuous Model Helps in Design of Functional Circuits
    '''
    def __init__(self, graph):
        super().__init__(graph)
        self.dxdt = self._squad(graph, gamma, h)








class SquadODE(ODE):

    def __init__(self, graph, transform, gamma=1, h=10):
        '''
        Transform a Activation/Inhibition Boolean network
        into an ODE system via the SQUAD transform.

        gamma -> decay rate
        h -> gain(?)

        If an int or float = same value for all variables
        other wise a dict with a default key and value and key value pairs for all others.
        '''
        super().__init__(graph)

        self.dxdt = self._squad(graph, gamma, h)
        print("done. ")



    def _squad(self, graph, gamma, h):

        self.n = graph.n

        self.gamma = parameter_to_array(gamma, self.n, graph.index)
        self.h = parameter_to_array(h, self.n, graph.index)

        # matrices
        self.Act, self.Inh = graph.primes_to_matrices()

        # needed for computations
        col_ones = np.ones((self.n))
        self._A1 = self.Act.dot(col_ones)
        self._a1 = (1+self._A1)/self._A1
        self._B1 = self.Inh.dot(col_ones)
        self._b1 = (1+self._B1)/self._B1

        dxdt = self._get_system()

        return dxdt

    def _omega(self, x):
        '''
        Equation 2 of Di Cara et al (2007)
        http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-462
        Based on Andre Blejec
        '''

        x = ensure_ndarray(x)

        col_ones = np.ones((self.n))

        Ax = self.Act.dot(x)
        a = ensure_ndarray( self._a1 * Ax/(1+Ax) )
        a[~np.isfinite(a)] = 1

        Bx = self.Inh.dot(x)
        b = ensure_ndarray( self._b1 * Bx/(1+Bx) )
        b[~np.isfinite(b)] = 0

        o = ensure_ndarray(a * (1-b))
        o[np.where(self._A1 + self._B1 ==0)] = 0

        return o

    def _dxdt_transform(self, x, w):
        return (  -np.exp(0.5*self.h) + np.exp(-self.h*(w-0.5))  ) / \
               (  (1-np.exp(0.5*self.h))*(1+np.exp(-self.h*(w-0.5)))  ) - self.gamma*x

    def _get_system(self, off_nodes=[]):

        def dxdt(t, x_array, *args):
            w = self._omega(x_array)
            d = self._dxdt_transform(x_array, w)
            for i in off_nodes:
                d[i] = 0
            return d

        return dxdt




ode_classes = {'squad':SquadODE,
               'shao':ShaoODE,
               'hillcube':BoolCube,
               'normalisedhillcube':BoolCube,
               'boolcube':BoolCube}