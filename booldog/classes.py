'''Plain dataclasses used throughout `booldog`: `BoolDogNode` (a single node's
identifier/name/rule) and `BoolDogModelInfo` (model metadata). A
`booldog.network.BoolDogModel` network is fundamentally a dict of
`BoolDogNode` objects, keyed by identifier.
'''

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

@dataclass
class BoolDogNode():
    '''A dataclass to represent a node in a Boolean network.'''
    identifier: str
    '''The identifier of the node.'''
    rule: str
    '''The Boolean rule of the node in Bnet format.'''

    # optional attributes
    name: Optional[str] = None
    '''The (optional) name of the node. Defaults to `identifier` if not given.'''

    def __post_init__(self):
        '''Default `name` to `identifier` if `name` was not given (falsy).'''
        if not self.name:
            self.name = self.identifier

    def __repr__(self):
        '''String representation showing the node's name and rule.'''
        return f"BoolDogNode(name='{self.name}', rule='{self.rule}')"


@dataclass
class BoolDogModelInfo():
    '''A dataclass to represent metadata about a Boolean network model.'''

    identifier: Optional[str] = None
    '''The identifier of the model. Defaults to "BoolDogModel" if not given
    (falsy). Set to the model's own id by some readers (e.g. the SBML-qual
    reader sets this to the SBML model id); used e.g. to title Cytoscape
    exports.'''
    source: Optional[str] = None
    '''The source of the model (e.g., file name, database, etc.).'''
    source_format: Optional[str] = None
    '''The exchange format the model was read from (e.g. "bnet", "primes",
    "sbml-qual"), as recorded by the corresponding `booldog.io` reader.'''
    notes: Optional[str] = None
    '''Free-text notes about the model. Not currently set or read anywhere
    else in the package.'''

    def __post_init__(self):
        '''Default `identifier` to "BoolDogModel" if not given (falsy).'''
        if not self.identifier:
            self.identifier = "BoolDogModel"

    def __repr__(self):
        '''String representation showing the model's identifier and source.'''
        return f"BoolDogModelInfo(identifier='{self.identifier}', source='{self.source}')"
