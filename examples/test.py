import sys
sys.path.insert(0,'/home/cbleker/research/NIB/squad/BoolDoG')
import booldog

from PyBoolNet import InteractionGraphs


interactions = booldog.utils.io.import_graphml('./networks/zagorscak-potatostress.graphml', 
                              inhibitor_symbol='t_shape')
                              
                              
                              
g = booldog.utils.io.SquadInteractions(interactions, default=1)                              
funcs = g.squad_update_funcs(interactions, default=1)


print("---------")

VPg=1
HCPro=0
CI=0
P3=0
PR11B25=0
PCD=0
bla_state = [VPg, HCPro, CI, P3, PR11B25, PCD]
print(funcs['virus'](*bla_state))

print("---------")

insect=1
CLH12JR1PR13THIVSP12=0
bla_state = [insect, CLH12JR1PR13THIVSP12]
print(funcs['insect'](*bla_state))

print("---------")
B = booldog.RegulatoryNetwork('./networks/zagorscak-potatostress.graphml', 
                              data_format='graphml', 
                              inhibitor_symbol='t_shape'
                             )

from PyBoolNet import InteractionGraphs
intgraph = InteractionGraphs.primes2igraph(B.primes)
print("num interactions ", intgraph.number_of_edges())
