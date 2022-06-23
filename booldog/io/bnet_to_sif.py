import re
import numpy as np


def write_sif

def clean_line(line):
    line = line.strip()
    if line == "targets, factors":
        return None

    m = re.match(r'^([^#]*)', line)
    if not m:
        return None

    return m.groups()[0].strip()


def parse_bnet(bnet):
    targets = set()

    with open(bnet, 'r') as f, open("test.tsv", "w") as out:
        for line in f:
            line = clean_line(line)
            if line:
                target, factors = line.split(",")
                target = target.strip()
                if target in targets:
                    print(f"{target} already has an update function")
                    continue

                targets.update([target])

                final, edges = resolve_brackets(factors)
                edges += [f"{final}\t{target}"]

                out.write('\n'.join(edges))
                out.write('\n')


    with open("test-annot.tsv", "w") as out:
        out.write("name\ttype\tlabel\n")
        for and_i in ands:
            out.write(f"{and_i}\tlogical and\tand\n")
        for or_i in ors:
            out.write(f"{or_i}\tlogical or\tor\n")
        for not_i in nots:
            out.write(f"{not_i}\tlogical not\tnot\n")

def resolve(s):
    '''
    parse a string without brackets

    Returns
    -------

    r: str
        Replacement string

    edges: list of str
        Edges from string
    '''

    global and_i, or_i, not_i

    new_s = []
    edges = []

    # if len(factors) == 1:
    #     return s, []

    for factor in s.split("|"):
        nodes = [s.strip() for s in factor.split("&")]
        for i, node in enumerate(nodes):
            m =  re.match(r".*?!(.*)", node)
            if m:
                node = m.groups()[0].strip()

                not_node = f"not_{not_i}"
                not_i += 1

                edges += [f"{node}\t{not_node}"]

                nodes[i] = not_node
                nots.append(not_node)

        if len(nodes) > 1:

            and_node = f"and_{and_i}"
            and_i += 1

            edges += [f"{node}\t{and_node}" for node in nodes]
            new_s.append(and_node)

            ands.append(and_node)

        else:
            new_s.append(nodes[0])

    if len(new_s) > 1:
        or_node = f"or_{or_i}"
        or_i += 1
        edges += [f"{node}\t{or_node}" for node in new_s]
        ors.append(or_node)


        return or_node, edges


    else:
        return new_s[0], edges


def resolve_brackets(s):

    groups = []
    levels = []

    group_index = 1
    num_open_brackets = 0
    group = '0'
    previous_groups = []
    for i, c in enumerate(s):

        if c == '(':
            num_open_brackets += 1
            previous_groups.append(group)
            group = f"{group}.{group_index}"

        elif c == ')':
            num_open_brackets += -1
            group = previous_groups.pop()

            group_index +=1

        groups.append(group)
        levels.append(num_open_brackets)

    groups = np.array(groups)
    levels = np.array(levels)

    print(groups)

    max_level = max(levels)

    groups_to_parse = np.unique([groups[i] for i in np.where(levels==max_level)[0]])

    print(groups_to_parse)


    starts = []
    ends = []
    new_sub_s = []

    new_s = ''

    all_edges = []
    for g in groups_to_parse:
        indices = np.nonzero(groups==g)[0]

        print(indices)

        if s[indices[0]] == "(":
            sub_s = ''.join(s[i] for i in indices[1:]).strip()
        else:
            sub_s = ''.join(s[i] for i in indices[:]).strip()
        replace, edges = resolve(sub_s)

        print(f'replace "{sub_s}" with "{replace}"')

        starts.append(indices[0])
        ends.append(indices[-1])
        new_sub_s.append(replace)

        all_edges += edges

    # print(starts, ends, new_sub_s)

    # ends.append(len(s))

    starts = np.array(starts)
    ends = np.array(ends)

    idx = np.argsort(starts)
    starts = starts[idx]
    ends = ends[idx]

    base_start = 0
    new_s = ''
    for i, (start, end) in enumerate(zip(starts, ends)):
        new_s += s[base_start:start]
        new_s += new_sub_s[i]
        base_start = end+2
    new_s  += s[base_start:]
    new_s = new_s.strip()

    print(new_s)
    if re.findall(r"[&|()]", new_s):
        new_s, new_edges = resolve_brackets(new_s)
        all_edges += new_edges


    return new_s, all_edges


# s = "( ADK.cyt.p & DZR.cyt.m | tZRMP.cyt.m | APT.cyt.p & DZ.cyt.m ) & ! (tZ.cyt.m & AHK234.cyt.p | !(cZ.cyt.m & AHK234.cyt.p) | DZ.cyt.m & !AHK234.cyt.p | iP.cyt.m & AHK234.cyt.p )"
# print(s)
# new_s, edges = resolve_brackets(s)

# print("----\nfinal", new_s)

# print(edges)

# with open("test.tsv", "w") as out:
#     out.write('\n'.join(edges))
#     out.write('\n')


# with open("test-annot.tsv", "w") as out:
#     out.write("name\ttype\tlabel\n")
#     for and_i in ands:
#         out.write(f"{and_i}\tlogical and\t&\n")
#     for or_i in ors:
#         out.write(f"{or_i}\tlogical or\tor\n")
#     for not_i in nots:
#         out.write(f"{not_i}\tlogical not\tnot\n")

# print(ors, ands, nots)



global and_i, or_i, not_i
and_i = 0
or_i = 0
not_i = 0

ands = []
ors = []
nots = []

import sys
bnet = sys.argv[1]

#bnet = "/home/cbleker/research/NIB/ADAPT/skm-webapp/downloads_data/pss-boolnet-restricted.bnet"
parse_bnet(bnet)