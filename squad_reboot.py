import numpy as np
import os
import sys

from collections import defaultdict
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# hide runtime warnings (divide by zero, multiply by inf)
# This is a bad idea...
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in multiply")

# https://github.com/hklarner/PyBoolNet/blob/master/Docs/Sphinx/source/Development.rst
PyBoolNet_path = os.path.join(os.path.abspath("."), "PyBoolNet/")
sys.path.insert(0, PyBoolNet_path)
from PyBoolNet import QuineMcCluskey as QMC
from PyBoolNet import AspSolver

class SquadRegulatoryNetwork:
    
    def __init__(self, g):
        
        self.boolean_graph = g # [here] add format translations?
        self.n = len(self.boolean_graph)
        
        # translate node string names to indices
        self.keys = {i:node for node, i in enumerate(self.boolean_graph.keys())} 
        
        # matrices
        self.Act, self.Inh = self._logic_matrices()

        # needed for computations
        col_ones = np.ones((self.n))
        self._A1 = self.Act.dot(col_ones)
        self._a1 = (1+self._A1)/self._A1
        self._B1 = self.Inh.dot(col_ones)
        self._b1 = (1+self._B1)/self._B1 
        
    def _ensure_ndarray(self, v):
        if not type(v) == np.ndarray:
            return np.array([*v])
        else:
            return v

    def _omega(self, x):
        '''Equation 2 of Di Cara et al (2007) 
           http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-462'''

        x = self._ensure_ndarray(x)

        col_ones = np.ones((self.n))

        Ax = self.Act.dot(x)
        a = self._ensure_ndarray( self._a1 * Ax/(1+Ax) )
        a[~np.isfinite(a)] = 1

        Bx = self.Inh.dot(x)
        b = self._ensure_ndarray( self._b1 * Bx/(1+Bx) )
        b[~np.isfinite(b)] = 0

        o = self._ensure_ndarray(a * (1-b))
        o[np.where(self._A1 + self._B1 ==0)] = 0 

        return o

    def _dxdt_transform(self, x, w, h, gamma):
        return (  -np.exp(0.5*h) + np.exp(-h*(w-0.5))  ) / \
               (  (1-np.exp(0.5*h))*(1+np.exp(-h*(w-0.5)))  ) - gamma*x

    def _get_system(self, h=50, gamma=1, off=[]):  

        def dxdt(t, x_array, *args):
            w = self._omega(x_array)
            d = self._dxdt_transform(x_array, w, h, gamma)
            for i in off:
                d[i] = 0
            return d

        return dxdt   

    def _event(self, t, x, event_t, *args):
        return t - event_t
    _event.terminal = True

    def _logic_matrices(self):
        # make up logic matrices
        # these can be made sparse matrices without much effort
        # see https://docs.scipy.org/doc/scipy/reference/sparse.html to use e.g. dok_matrix 
        Act = np.zeros((self.n, self.n)) # activators
        Inh = np.zeros((self.n, self.n)) # inhibitors
        for node, d in self.boolean_graph.items():
            for other_node, sign in d.items():
                if sign == "+":
                    Act[self.keys[node], self.keys[other_node]] = 1
                elif sign == "-":
                    Inh[self.keys[node], self.keys[other_node]] = 1
                else:
                    print("Warning: Issue with edge: ", node, other_node)
        return self._ensure_ndarray(Act), self._ensure_ndarray(Inh)

    def _get_node_bool_func(self, args, node):
        def func(*func_input):        
            state = np.ones(n)
            for other_node, other_node_state in zip(args, func_input):
                state[keys[other_node]] = other_node_state

            col_ones = np.ones((n))
            if Inh[keys[node],:].dot(col_ones):
                inh = Inh[keys[node],:].dot(state) > 0
            else:
                inh = 0
            if Act[keys[node],:].dot(col_ones):
                act = Act[keys[node],:].dot(state) > 0
            else:
                act = 1
            node_state = act * (1-inh)

            #print(node, state, inh, act, node_state)
            return node_state

        func.depends = [pyboolvar[node] for node in args]
        
        return func

    def _rename_variables(self, orig_d, rev_d, direction):
        
        if direction == 'forward':
            pyboolvar = {node:node + "var" for node in orig_d.keys()}
            pyboolvar_rev = {node + "var":node for node in orig_d.keys()}
            return pyboolvar, pyboolvar_rev
        
        elif direction == 'backward':
            orig_d = {}
            for node, d in d.items():
                node = rev_d[node]
                p_rev = []
                for sublist in d:
                    sublist_rev = []
                    for sub_d in sublist:
                        sub_d_rev = {}
                        for other_node, state in sub_d.items():
                            other_node = rev_d[other_node]
                            sub_d_rev[other_node] = state
                        sublist_rev.append(sub_d_rev)
                    p_rev.append(sublist_rev)
                orig_d[node] = p_rev
            return orig_d

    def bool_steady_states(self):
        
        # need to change variable names because pf pyboolnet errors
        tmp_keys, tmp_key_rev = self._rename_variables(self.keys(), {}, 'forward')

        funcs = {tmp_keys[node]:1 for node in self.boolean_graph.keys()}
        for node, d in self.boolean_graph.items():
            print(node)
            args = list(d.keys())
            func = self._get_node_bool_func(args, node)
            funcs[tmp_keys[node]] = func
        primes = QMC.functions2primes(funcs)
        primes = self._rename_variables(primes, tmp_key_rev, 'backward')
        
        steady_states = AspSolver.steady_states(primes)
        return steady_states
        

    def dynamic_simulation(self, 
                           events={}, 
                           t_min=0, t_max=30, 
                           initial_state=0, 
                           plot=True,
                           gamma=1, h=50):
        
        all_results = []

        initial_state_array = np.zeros(self.n)
        if not (initial_state == 0):
            for node, value in initial_state.items():
                initial_state_array[self.keys[node]] = value


        # fetch the system of eqn
        dxdt = self._get_system(h=h, gamma=gamma)

        # 1. copy events to new dict
        # 2. set a duration for all events
        # 3. add end of each pertubation as an event (so that dx/dt of the node can be reset)
        events_d = defaultdict(dict)
        for event_t  in events.keys():
            for node, perturb in events[event_t].items():
                events_d[event_t][node] = perturb
                events_d[event_t][node]["reset"] = False
                if not ("duration" in perturb):
                    perturb_duration = events_d[event_t][node]["duration"] = 0
                else:
                    perturb_duration = perturb["duration"]

                if perturb_duration > 0:                
                    perturb_end = min(t_max, int(event_t + perturb_duration))
                    events_d[perturb_end] = {node:{"reset":True}}

        # add a last event to get to end of simulation
        # i.o.t. avoid using `while True`
        # # https://stackoverflow.com/a/50703835/4996681
        events_d[t_max] = {} 

        off_nodes = set()
        for event_t  in sorted(events_d.keys()):
            result = solve_ivp(fun=dxdt, t_span=(t_min, t_max), y0=initial_state_array, 
                                         events=self._event, args=(event_t, self.Act.T, self.Inh.T, h, gamma), 
                                         dense_output=True, max_step=0.01)
            all_results.append(result)
            current_t = result.t[-1]

            if (result.status == 0) or (current_t == t_max):
                print("Status: End")
                break
            else:
                s = "Status: Event at {t:3d}:\n".format(t=int(current_t))

                t_min = current_t
                initial_state_array = result.y[:, -1].copy()
                for node, perturb in events_d[event_t].items():

                    if perturb["reset"]:
                        # ending a pertubation (put back dx/dt of this node)
                        s += " "*10 + "{node} -> released\n".format(node=node)
                        off_nodes.remove(self.keys[node])

                    else:
                        # starting a pertubation
                        initial_state_array[self.keys[node]] = perturb['perturbation']
                        s += " "*10 + "{node} -> {perturb_value:.2f} (duration {duration:2d})\n".format(\
                                                           node=node,                              
                                                           perturb_value=perturb['perturbation'], duration=perturb["duration"])

                        if perturb['duration'] > 0:
                            # dt/dt of this node should be 0 for "duration"
                            off_nodes.add(self.keys[node])
                            dxdt = self._get_system(h=h, gamma=gamma, off=off_nodes)

                print(s)



        combined_t = np.concatenate([result.t for result in all_results])
        combined_y = np.concatenate([result.y for result in all_results], axis=1).T

        if plot: 
            lines = plt.plot(combined_t, combined_y, '-')
            plt.legend(lines, self.boolean_graph.keys())