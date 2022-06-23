''' Read Boolean networks from various formats and return a RegulatoryNetwork
object.
'''

from .read import read_primes, read_bnet, read_interactions, read_graphml





# sphinx stuff
__all_exports = [read_primes, read_bnet, read_interactions, read_graphml]

for e in __all_exports:
    e.__module__ = __name__

__all__ = [e.__name__ for e in __all_exports]



# __all__ = ['read_primes', 'read_bnet', 'read_interactions', 'read_graphml', 'read']