'''Examples for testing.'''

class BooleanNetworkExamples:

    BNET = """targets, factors
    node_A, (node_B & !node_E) | (!node_B & node_E)
    node_C, node_A & node_B
    node_E, !node_C | node_D
    node_D, node_D
    """

    PRIMES = {
        'node_A': [[{'node_B': 0, 'node_E': 0}, {'node_B': 1, 'node_E': 1}], [{'node_B': 0, 'node_E': 1}, {'node_B': 1, 'node_E': 0}]],
        'node_B': [[{'node_B': 0}], [{'node_B': 1}]],
        'node_C': [[{'node_B': 0}, {'node_A': 0}], [{'node_A': 1, 'node_B': 1}]],
        'node_D': [[{'node_D': 0}], [{'node_D': 1}]],
        'node_E': [[{'node_C': 1, 'node_D': 0}], [{'node_D': 1}, {'node_C': 0}]]
    }

    SBMLQUAL_FILE = "data/example.xml"

    TABULARQUAL_FILE = "data/example.xlsx"

class InteractionNetworkExamples:

    PRIMES_SQUAD = {
        'node_A': [[{'node_E': 1}, {'node_B': 1}], [{'node_B': 0, 'node_E': 0}]],
        'node_B': [[{'node_B': 0}], [{'node_B': 1}]],
        'node_C': [[{'node_A': 0, 'node_B': 0}], [{'node_B': 1}, {'node_A': 1}]],
        'node_D': [[{'node_D': 0}], [{'node_D': 1 }]],
        'node_E': [[{'node_D': 0}, {'node_C': 1 }], [{'node_C': 0, 'node_D': 1 }]]}

    INTERACTIONS = [
        ('node_A', 'node_C', '+'),
        ('node_B', 'node_A', '-'),
        ('node_B', 'node_B', '+'),
        ('node_B', 'node_C', '+'),
        ('node_C', 'node_E', '-'),
        ('node_D', 'node_D', '+'),
        ('node_D', 'node_E', '+'),
        ('node_E', 'node_A', '-')
    ]

    DICT_OF_DICT = {
        'node_A': {'node_C': {'interaction': '+'}},
        'node_B': {'node_A': {'interaction': '-'},        'node_B': {'interaction': '+'},        'node_C': {'interaction': '+'}},
        'node_C': {'node_E': {'interaction': '-'}},
        'node_D': {'node_D': {'interaction': '+'},        'node_E': {'interaction': '+'}},
        'node_E': {'node_A': {'interaction': '-'}}
    }

    SIF_FILE = "data/example.sif"

    GRAPHML_FILE = "data/example.graphml"

    GRAPHML_YED_FILE = "data/example.yEd.graphml"
