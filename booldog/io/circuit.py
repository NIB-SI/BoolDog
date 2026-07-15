'''
Convert a BoolDogModel's Boolean rules into a logic-circuit graph.

Each node's Boolean rule (a string using ``&``, ``|``, ``!``) is parsed and
expanded into explicit "logical" AND/OR/NOT nodes wired up with edges, in
addition to the original "species" nodes. This produces a
:class:`networkx.DiGraph` circuit representation of the model, as used by
``booldog2networkx``/``booldog2igraph`` when ``as_logic_circuit=True``.
'''
import logging

import re
import numpy as np
import networkx as nx
from pyboolnet.interaction_graphs import primes2igraph

logger = logging.getLogger(__name__)

##############################
# Utility functions
##############################

PAREN_RE = re.compile(r"\([^()]*\)")

def clean_line(line):
    ''' Strip a bnet-style line down to its Boolean-rule content.

    Parameters
    ----------
    line : str
        A single line of text (e.g. from a bnet file).

    Returns
    -------
    str or None
        ``None`` if, once stripped, `line` is exactly the bnet header line
        ("targets, factors"). Otherwise, the part of `line` before any
        trailing "#" comment, stripped of surrounding whitespace.
    '''
    line = line.strip()
    if line == "targets, factors":
        return None

    m = re.match(r'^([^#]*)', line)
    if not m:
        return None

    return m.groups()[0].strip()


class BooleanOperators:
    ''' A simple class to define the Boolean operators as class attributes. This
    is used to make the code more readable and maintainable.
    '''
    AND = "&"
    '''str : symbol used for the Boolean AND operator in rule strings.'''
    OR = "|"
    '''str : symbol used for the Boolean OR operator in rule strings.'''
    NOT = "!"
    '''str : symbol used for the Boolean NOT operator in rule strings.'''

class BooleanDiGraph(nx.DiGraph):
    ''' A :class:`networkx.DiGraph` subclass that also tracks how many AND, OR,
    and NOT logic nodes it contains, so that newly added logic nodes get
    unique, sequential labels (e.g. "and_0", "and_1", ...).
    '''

    def __init__(self):

        super().__init__()

        self.and_count = 0
        '''int : number of AND nodes added so far, used to generate the next unique "and_<n>" node id.'''
        self.or_count = 0
        '''int : number of OR nodes added so far, used to generate the next unique "or_<n>" node id.'''
        self.not_count = 0
        '''int : number of NOT nodes added so far, used to generate the next unique "not_<n>" node id.'''

    def add_and(self):
        '''
        Add a new AND logic node ("and_<n>") to the graph.

        Returns
        -------
        node : str
            The identifier of the newly added node (e.g. "and_0").
        '''

        node = f"and_{self.and_count}"
        self.add_node(node, type="logical_and", label="and")
        self.and_count += 1

        return node

    def add_or(self):
        '''
        Add a new OR logic node ("or_<n>") to the graph.

        Returns
        -------
        node : str
            The identifier of the newly added node (e.g. "or_0").
        '''

        node = f"or_{self.or_count}"
        self.add_node(node, type="logical_or", label="or")
        self.or_count += 1

        return node

    def add_not(self):
        '''
        Add a new NOT logic node ("not_<n>") to the graph.

        Returns
        -------
        node : str
            The identifier of the newly added node (e.g. "not_0").
        '''

        node = f"not_{self.not_count}"
        self.add_node(node, type="logical_not", label="not")
        self.not_count += 1

        return node

##############################
# Exchange functions
##############################


def booldog2circuit(model):
    '''Export a BoolDog Boolean model to Networkx DiGraph

    Parameters
    ----------
    model : booldog:BoolDogModel
        A BoolDog object representing a Boolean network.

    Returns
    -------
    graph : networkx.DiGraph
        A plain networkx.DiGraph (built from an internal BooleanDiGraph) with
        one "species" node per model node, plus additional "logical" AND/OR/NOT
        nodes and edges representing each node's Boolean rule. An edge of type
        "update_function" connects the final logic node of a rule to the
        species node it updates. Nodes whose rule is empty/falsy have no
        incoming update-function edge (a warning is logged instead).
    '''


    graph = BooleanDiGraph()
    for node_id in model.nodes:
        graph.add_node(node_id, type="species")

    def update_digraph(nodes, operator):

        if operator == BooleanOperators.NOT:
            if len(nodes) != 1:
                raise ValueError("NOT operator can only be applied to a single node.")

            not_node = graph.add_not()
            graph.add_edge(nodes[0], not_node, type="to_not")
            return not_node

        if operator == BooleanOperators.AND:
            and_node = graph.add_and()
            for node in nodes:
                graph.add_edge(node, and_node, type="to_and")
            return and_node

        if operator == BooleanOperators.OR:
            or_node = graph.add_or()
            for node in nodes:
                graph.add_edge(node, or_node, type="to_or")
            return or_node

    for target_id in model.nodes:
        target = model.nodes[target_id]
        rule = target.rule

        if not rule:
            logger.warning("No update function for %s.", target)
            continue

        logger.debug("Processing rule for %s: %s", target.identifier, rule)
        final_node = resolve_rule(rule, update_digraph)
        graph.add_edge(final_node, target.identifier, type="update_function")

    return nx.DiGraph(graph)

def resolve_rule(rule, processor):
    '''
    Parse a Boolean rule and replace it with "logical" nodes and edges in the
    graph.

    The rule can be nested with brackets, which are resolved first.

    Parameters
    ----------
    rule : str
        A Boolean rule (e.g. "A & B | C")
    processor: function
        A function that takes a list of nodes and an operator, processes them as
        needed. The operator can be "AND", "OR", or "NOT". The function should
        return the new node representing the result of applying the operator to the
        nodes.

    Returns
    -------
    new_s : str
        Replacement string (the final node)

    '''

    rule = rule.strip()

    # Get the positions of the innermost parentheses in a string.
    while True:
        match = PAREN_RE.search(rule)
        if not match:
            break

        sub_expr = match.group(0)[1:-1].strip() # remove parentheses
        replacement = resolve_factor(sub_expr, processor)
        logger.debug("Replace %s with %s", sub_expr, replacement)

        rule = rule[:match.start()] + replacement + rule[match.end():] # replace parentheses

    # remaining (no parentheses)
    if re.search(r"[&|!]", rule):
        rule = resolve_factor(rule, processor)

    return rule

def resolve_factor(factor, processor):
    '''
    Parse a single factor of a Boolean rule.

    The factor can be a list of nodes connected by '&', '|', and '!' operators.
    The factor does not contain brackets.

    Parameters
    ----------
    factor : str
        A Boolean factor e.g. "A & B | C"

    processor: function
        A function that takes a list of nodes and an operator, processes them as
        needed. The operator can be "AND", "OR", or "NOT". The function should
        return the new node representing the result of applying the operator to the
        nodes.

    Returns
    -------
    s : str
        Replacement string (the final node representing this factor)

    '''

    new_s = []

    for term in factor.split("|"):
        nodes = [s.strip() for s in term.split("&")]

        # replace nodes with nots
        for i, node in enumerate(nodes):
            m = re.match(r".*?!(.*)", node)
            if m:
                node = m.groups()[0].strip()
                not_node = processor([node], BooleanOperators.NOT)
                nodes[i] = not_node

        if len(nodes) > 1:
            and_node = processor(nodes, BooleanOperators.AND)
            new_s.append(and_node)

        else:
            new_s.append(nodes[0])

    if len(new_s) > 1:
        return processor(new_s, BooleanOperators.OR)
    return new_s[0]
