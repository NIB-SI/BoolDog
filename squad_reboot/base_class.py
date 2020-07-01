import numpy as np
import os
import sys

import matplotlib.pyplot as plt

from collections import defaultdict
from scipy.integrate import solve_ivp

# hide runtime warnings (divide by zero, multiply by inf)
# This is a bad idea...
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
warnings.filterwarnings("ignore", message="invalid value encountered in multiply")

from PyBoolNet import QuineMcCluskey as QMC
from PyBoolNet import AspSolver


from .utils import *

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
        b = _ensure_ndarray( self._b1 * Bx/(1+Bx) )
        b[~np.isfinite(b)] = 0

        o = ensure_ndarray(a * (1-b))
        o[np.where(self._A1 + self._B1 ==0)] = 0 

        return o

    def _dxdt_transform(self, x, w, h, gamma):
        return (  -np.exp(0.5*h) + np.exp(-h*(w-0.5))  ) / \
               (  (1-np.exp(0.5*h))*(1+np.exp(-h*(w-0.5)))  ) - gamma*x
   
    def _get_system(self, h, gamma, off=[]):  

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
        return ensure_ndarray(Act), ensure_ndarray(Inh)

    def _get_node_bool_func(self, args, node, tmp_keys):
        def func(*func_input):        
            state = np.ones(self.n)
            for other_node, other_node_state in zip(args, func_input):
                state[self.keys[other_node]] = other_node_state

            col_ones = np.ones((self.n))
            if self.Inh[self.keys[node],:].dot(col_ones):
                inh = self.Inh[self.keys[node],:].dot(state) > 0
            else:
                inh = 0
            if self.Act[self.keys[node],:].dot(col_ones):
                act = self.Act[self.keys[node],:].dot(state) > 0
            else:
                act = 1
            node_state = act * (1-inh)

            #print(node, state, inh, act, node_state)
            return node_state

        func.depends = [tmp_keys[node] for node in args]
        
        return func

    def _rename_variables(self, orig_d, rev_d, direction):
        
        if direction == 'forward':
            pyboolvar = {node:node + "var" for node in orig_d.keys()}
            pyboolvar_rev = {node + "var":node for node in orig_d.keys()}
            return pyboolvar, pyboolvar_rev
        
        elif direction == 'backward':
            new_d = {}
            for node, d in orig_d.items():
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
                new_d[node] = p_rev
            return new_d

    def bool_steady_states(self):
        
        # need to change variable names because of pyboolnet errors
        tmp_keys, tmp_key_rev = self._rename_variables(self.keys, {}, 'forward')

        funcs = {tmp_keys[node]:1 for node in self.boolean_graph.keys()}
        for node, d in self.boolean_graph.items():
            args = list(d.keys())
            func = self._get_node_bool_func(args, node, tmp_keys)
            funcs[tmp_keys[node]] = func
        primes = QMC.functions2primes(funcs)
        primes = self._rename_variables(primes, tmp_key_rev, 'backward')
        
        steady_states = AspSolver.steady_states(primes)
        return steady_states
        

    def dynamic_simulation(self, 
                           node_events={}, 
                           edge_events={},
                           t_min=0, t_max=30, 
                           initial_state=0, 
                           plot=True,
                           gamma=1, h=10):
        '''
        gamma -> decay rate
        h -> gain(?)
        '''
        
        all_results = []

        initial_state_array = np.zeros(self.n)
        if not (initial_state == 0):
            for node, value in initial_state.items():
                initial_state_array[self.keys[node]] = value
        
        gamma_array = np.ones(self.n)
        if type(gamma) in [int, float]:
            gamma_array = gamma_array * gamma
        else:
            if 'default' in gamma.keys():
                gamma_array = gamma_array * gamma['default']
            for node, value in gamma.items():
                if node == 'default': 
                    pass
                else:
                    gamma_array[self.keys[node]] = value            
        
        h_array = np.ones(self.n)
        if type(h) in [int, float]:
            h_array = h_array * h
        else: 
            if 'default' in h.keys():
                h_array = h_array * h['default']
            for node, value in h.items():
                if node == 'default': 
                    pass
                else:                
                    h_array[self.keys[node]] = value                   
                
        # fetch the system of eqn
        dxdt = self._get_system(h=h_array, gamma=gamma_array)

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
                        dxdt = self._get_system(h=h_array, gamma=gamma_array, off=off_nodes)

                    else:
                        # starting a pertubation
                        initial_state_array[self.keys[node]] = perturb['perturbation']
                        s += " "*10 + "{node} -> {perturb_value:.2f} (duration {duration:2d})\n".format(\
                                                           node=node,                              
                                                           perturb_value=perturb['perturbation'], 
                                                           duration=perturb["duration"])

                        if perturb['duration'] > 0:
                            # dt/dt of this node should be 0 for "duration"
                            off_nodes.add(self.keys[node])
                            dxdt = self._get_system(h=h_array, gamma=gamma_array, off=off_nodes)
                print(s)

        combined_t = np.concatenate([result.t for result in all_results])
        combined_y = np.concatenate([result.y for result in all_results], axis=1).T

        if plot: 
            lines = plt.plot(combined_t, combined_y, '-')
            plt.legend(lines, self.boolean_graph.keys())
            
        return combined_t, combined_y