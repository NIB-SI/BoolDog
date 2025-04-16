'''
Additional functions to parse bnet format

'''
import logging

import re
import numpy as np
import networkx as nx

logger = logging.getLogger(__name__)

##############################
# Utility functions
##############################


def clean_line(line):
    ''' Makes sure line is a Boolean rule, and not a header or comment
    '''
    line = line.strip()
    if line == "targets, factors":
        return None

    m = re.match(r'^([^#]*)', line)
    if not m:
        return None

    return m.groups()[0].strip()


class NetworkxBnet(nx.DiGraph):
    ''' Simple helper class to track logic nodes '''

    def __init__(self):

        super().__init__()

        self.and_count = 0
        self.or_count = 0
        self.not_count = 0

    def add_and(self):

        node = f"and_{self.and_count}"
        self.add_node(node, type="logical_and", label="and")
        self.and_count += 1

        return node

    def add_or(self):

        node = f"or_{self.or_count}"
        self.add_node(node, type="logical_or", label="or")
        self.or_count += 1

        return node

    def add_not(self):

        node = f"not_{self.not_count}"
        self.add_node(node, type="logical_not", label="not")
        self.not_count += 1

        return node


def booldog2networkx(network):
    '''
    Loads a graph from bnet into a networkx DiGraph
    '''

    graph = NetworkxBnet()

    for node in sorted(network.nodes):
        graph.add_node(node, type="species", label=node)

    for line in network.to_bnet(header=False).split("\n"):
        line = clean_line(line)
        if line:
            logger.debug("Working on line: %s", line)
            target, rule = target, rule = [
                s.strip() for s in line.split(",", 1)
            ]
            if rule:
                final_node = resolve_rule(rule, graph)
                graph.add_edge(final_node, target)

    return nx.DiGraph(graph)


def resolve_rule(rule, graph):

    groups = []
    levels = []

    group_index = 1
    num_open_brackets = 0
    group = '0'
    previous_groups = []
    for i, c in enumerate(rule):

        if c == '(':
            num_open_brackets += 1
            previous_groups.append(group)
            group = f"{group}.{group_index}"

        elif c == ')':
            num_open_brackets += -1
            group = previous_groups.pop()

            group_index += 1

        groups.append(group)
        levels.append(num_open_brackets)

    groups = np.array(groups)
    levels = np.array(levels)

    max_level = max(levels)

    groups_to_parse = np.unique(
        [groups[i] for i in np.where(levels == max_level)[0]])

    starts = []
    ends = []
    new_sub_s = []

    new_s = ''

    for g in groups_to_parse:
        indices = np.nonzero(groups == g)[0]

        if rule[indices[0]] == "(":
            sub_s = ''.join(rule[i] for i in indices[1:]).strip()
        else:
            sub_s = ''.join(rule[i] for i in indices[:]).strip()
        replace = resolve_factor(sub_s, graph)

        logger.debug("Replace %s with %s", sub_s, replace)

        starts.append(indices[0])
        ends.append(indices[-1])
        new_sub_s.append(replace)

    starts = np.array(starts)
    ends = np.array(ends)

    idx = np.argsort(starts)
    starts = starts[idx]
    ends = ends[idx]

    base_start = 0
    new_s = ''
    for i, (start, end) in enumerate(zip(starts, ends)):
        new_s += rule[base_start:start]
        new_s += new_sub_s[i]
        base_start = end + 2
    new_s += rule[base_start:]
    new_s = new_s.strip()

    logger.debug(new_s)
    if re.findall(r"[&|()]", new_s):
        new_s = resolve_rule(new_s, graph)

    return new_s

def resolve_factor(factor, graph):
    '''
    parse a string without brackets

    Returns
    -------

    r: str
        Replacement string (the final node)

    '''

    new_s = []

    for term in factor.split("|"):
        nodes = [s.strip() for s in term.split("&")]
        for i, node in enumerate(nodes):
            m = re.match(r".*?!(.*)", node)
            if m:
                node = m.groups()[0].strip()

                not_node = graph.add_not()
                graph.add_edge(node, not_node)

                # replace node with not
                nodes[i] = not_node

        if len(nodes) > 1:
            and_node = graph.add_and()
            for node in nodes:
                graph.add_edge(node, and_node)
            new_s.append(and_node)

        else:
            new_s.append(nodes[0])

    if len(new_s) > 1:
        or_node = graph.add_or()
        for node in new_s:
            graph.add_edge(node, or_node)
        return or_node

    return new_s[0]
