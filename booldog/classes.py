from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

@dataclass
class BoolDogNode():
    '''A dataclass to represent a node in a Boolean network.

    Attributes
    ----------
    identifier : str
        The identifier of the node.

    name : str
        The name of the node.

    rule : str
        The Boolean rule of the node in Bnet format.

    '''
    identifier: str
    rule: str

    # optional attributes
    name: Optional[str] = None

    def __post_init__(self):
        if not self.name:
            self.name = self.identifier

    def __repr__(self):
        return f"BoolDogNode(name='{self.name}', rule='{self.rule}')"


# class BoolDogRule():
#     '''A class to represent a Boolean rule in a Boolean network.

#     Attributes
#     ----------
#     target : str
#         The target node of the rule.

#     factors : list of str
#         The list of factor nodes in the rule.

#     operator : str
#         The logical operator of the rule (e.g., "AND", "OR", "NOT").

#     '''
#     def __init__(self, target, factors, operator):
#         self.target = target
#         self.factors = factors
#         self.operator = operator

#     def __repr__(self):
#         return f"BoolDogRule(target='{self.target}', operator='{self.operator}')"

@dataclass
class BoolDogModelInfo():
    '''A dataclass to represent metadata about a Boolean network model.

    Attributes
    ----------
    source : str
        The source of the model (e.g., file name, database, etc.).
    '''

    identifier: Optional[str] = None
    source: Optional[str] = None
    source_format: Optional[str] = None
    notes: Optional[str] = None

    def __post_init__(self):
        if not self.identifier:
            self.identifier = "BoolDogModel"

    def __repr__(self):
        return f"BoolDogModelInfo(identifier='{self.identifier}', source='{self.source}')"