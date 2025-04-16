'''Read-write SBML-qual files

Notes
=====

SBML-qual read code is inspired by:
* loadSBML.R of BoolNet (January 2024)
    https://github.com/cran/BoolNet/blob/af3a714c5bfa72ee7507db9c4eaf90ba2cd91809/R/loadSBML.R
* sbmlqual.py of CellNOpt (January 2024)
    https://github.com/cellnopt/cellnopt/blob/319c956e890ed7ae88c8e6e540de7061acfaae00/cno/io/sbmlqual.py

The function `parse_mathml` (and its dependencies) was initially written by
ChatGPT, and then fixed by hand.

'''
import re
import logging

try:
    import libsbml
    from libsbml import SBMLReader, SBMLNamespaces, SBMLDocument
    _SBML_AVAILABLE = True
except ImportError as e:
    _SBML_AVAILABLE = False

from booldog.utils import functions2mindnf

logger = logging.getLogger(__name__)

# libsbml sets certain missing values to SBML_MAX_INT, so to determine if they
# are missing, need to check against SBML_MAX_INT. Value from libsbml C code.
# TODO: is there a better option?
SBML_INT_MAX = 2147483647


class BoolDogSBMLException(Exception):
    '''Custom Exception for SBML parsing '''


def sbmlqual2bnet(file):
    '''Read SBML-qual xml file to a BooleanNetwork object, via bnet.

    Parameters
    ----------
    file : str or path-like
        Path to SBML-qual file containing a Boolean network.

    Returns
    -------
    bnet: str
        bnet representation of the Boolean network.
    '''

    sbmldoc = SBMLReader().readSBML(file)

    # Check for errors
    num_errors = sbmldoc.getNumErrors()
    if num_errors > 0:
        print(f"SBML file has {num_errors} error(s):")
        for i in range(num_errors):
            print(f"  {sbmldoc.getError(i).getMessage()}")
        raise ValueError("SBML file has errors, not loaded.")

    model = sbmldoc.getModel()

    # Extract 'qual' plugin
    have_plugin = False
    for i in range(sbmldoc.getNumPlugins()):
        p = model.getPlugin(i)
        if p.getPackageName() == 'qual':
            have_plugin = True
            break
    if not have_plugin:
        raise ValueError("SBML file does not have 'qual' plugin, not loaded.")

    all_species = []
    for i in range(p.getNumQualitativeSpecies()):
        species = p.getQualitativeSpecies(i)
        all_species.append(species.id)

    bnet = "targets, factors\n"
    for i in range(p.getNumTransitions()):
        transition = p.getTransition(i)
        inputs, outputs = parse_transition_io(transition, all_species)
        logic_rule = parse_transition_function(transition, all_species, inputs)

        for o in outputs:
            bnet += f"{o['species']}, {logic_rule}\n"

    return bnet


def parse_transition_io(transition, all_species):
    '''Extract transition inputs and outputs.

    Parameters
    ----------
    transition : libsbml::Transition
        Transition specifying the logical rule associated with the
        Transition outputs.
    all_species: list
        List ids of all species present in model

    Returns
    -------
    inputs: dict
        Dictionary of id: species information for this transition's inputs
    outputs: list
        List of this transition's outputs (dependents)
    '''

    inputs = {}
    for sp in transition.getListOfInputs():
        d = {
            "id": sp.getId(),
            "species": sp.getQualitativeSpecies(),
            "sign": sp.getSign(),
            "threshold": None if sp.getThresholdLevel() == SBML_INT_MAX else sp.getThresholdLevel(),
            "transition_effect": sp.getTransitionEffect()
        }

        if not (d["species"] in all_species):
            logger.warning("Species '%s' not defined in model", d["species"])

        if d["transition_effect"] != libsbml.INPUT_TRANSITION_EFFECT_NONE:
            logger.warning("Effect '%s' is not 'None' for transition input %s",
                           d["transition_effect"], d["id"])

        if not (d["threshold"] in (0, 1)):
            logger.warning(
                "Threshold '%s' is not 0 or 1 for transition input %s",
                d["threshold"], d["id"])

        inputs[d["id"]] = d

    outputs = []
    for sp in transition.getListOfOutputs():
        d = {
            "species": sp.getQualitativeSpecies(),
            "transition_effect": sp.getTransitionEffect(),
            "output_level": None if sp.getOutputLevel() == SBML_INT_MAX else sp.getOutputLevel()
        }
        if not (d["species"] in all_species):
            logger.warning("Species '%s' not defined in model", d["species"])
            break

        if d["transition_effect"] != libsbml.OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL:
            logger.warning(
                "Transition effect '%s' not defined for Boolean model",
                d["transition_effect"])
            break

        if d["output_level"] is not None:
            logger.warning("Output level '%s' not supported for Boolean model",
                           d["output_level"])
            break

        outputs.append(d)

    return inputs, outputs


def parse_transition_function(transition, all_species, inputs):
    '''Parse transition function to a logical rule.

    Parameters
    ----------
    transition : libsbml::Transition
        Transition specifying the logical rule associated with the
        Transition outputs.
    all_species: list
        List ids of all species present in model
    inputs: dict
        Dictionary of id: species information for this transition's inputs

    Returns
    -------
    logic_rule: str
        Logic rule of this transition (in bnet format)
    '''

    function_terms = transition.getListOfFunctionTerms()

    activation_terms = []
    inhibition_terms = []

    for function in function_terms:

        output_level = function.getResultLevel()
        if not output_level in (0, 1):
            logger.warning('ResultLevel %i is not 0 or 1', output_level)

        if function.name.lower() == "defaultterm":
            continue

        mathml = function.getMath()

        try:
            logic_rule = parse_MathML(mathml, all_species, inputs)
        except Exception as err:
            raise BoolDogSBMLException(
                f'Fatal error occurred in parsing  transition {transition.getId()}'
            ) from err

        logger.debug('%s <==> %s', libsbml.formulaToL3String(mathml),
                     logic_rule)

        if output_level == 1:
            activation_terms.append(logic_rule)
        else:
            inhibition_terms.append(logic_rule)

    default_term = transition.getDefaultTerm().getResultLevel()
    # The DefaultTerm defines the default result of a Transition. This term
    # is used if there are no other FunctionTerm elements or if none of the
    # Math elements of the FunctionTerm elements evaluates to “true”.

    # BoolNet logic, not sure I agree...
    # if default_term == 0:
    #     # default is off, can only activate
    #     logic_rule = "( " + " | ".join(activation_terms) + " )"
    #     if inhibition_terms:
    #         logger.info("Ignoring contradictions in rule for %s", ','.join(inputs))
    # else:
    #     # default is on, can only inhibt
    #     logic_rule = "!( " + " | ".join(inhibition_terms) + " )"
    #     if activation_terms:
    #         logger.info("Ignoring contradictions in rule for %s", ','.join(inputs))

    # alternative logic...
    if activation_terms:
        logic_rule = "( " + " | ".join(activation_terms) + " )"
        if inhibition_terms:
            logic_rule += "& !( " + " | ".join(inhibition_terms) + " )"

    elif inhibition_terms:
        logic_rule = "!( " + " | ".join(inhibition_terms) + " )"

    else:
        logic_rule = default_term

    return logic_rule


def handle_xor(children):
    '''Create disjunctive normal form (DNF) of an xor '''

    def xor(*l):
        return not (sum(l) % 2 == 0)

    xor.depends = children
    mindnf = functions2mindnf({"xor_function": xor})['xor_function']

    return mindnf


def handle_not(children):
    '''Create disjunctive normal form (DNF) of an xor '''

    if len(children) != 1:
        raise ValueError("Unary operator \"not\" can only have one argument!")

    return f"! {children[0]}"


def handle_comparison(operator, children):
    '''Return bnet form of an (mathml) operator between two children.'''

    if len(children) != 2:
        raise ValueError(f"Operator \"{operator}\" requires two operands!")

    is_const = [isinstance(child, int) for child in children]

    if all(is_const):
        children = [int(child) for child in children]
        return str({
            "eq": lambda x, y: x == y,
            "neq": lambda x, y: x != y,
            "gt": lambda x, y: x > y,
            "lt": lambda x, y: x < y,
            "geq": lambda x, y: x >= y,
            "leq": lambda x, y: x <= y,
        }[operator](children[0], children[1]))

    if any(is_const):
        const_child = int(children[is_const.index(True)])
        var_child = children[is_const.index(False)]
        return {
            "eq": f"{var_child}" if const_child == 1 else f"!{var_child}",
            "neq": f"{var_child}" if const_child == 0 else f"!{var_child}",
            "gt": "0" if const_child == 1 else f"{var_child}",
            "lt": "0" if const_child == 0 else f"!{var_child}",
            "geq": "1" if const_child == 0 else f"{var_child}",
            "leq": "1" if const_child == 1 else f"!{var_child}",
        }[operator]

    return {
        "eq":
        f"(({children[0]} & {children[1]}) | (!{children[0]} & !{children[1]}))",
        "neq":
        f"(({children[0]} & !{children[1]}) | (!{children[0]} & {children[1]}))",
        "gt": f"({children[0]} & !{children[1]})",
        "lt": f"(!{children[0]} & {children[1]})",
        "geq": f"({children[0]} | !{children[1]})",
        "leq": f"(!{children[0]} | {children[1]})",
    }[operator]


def parse_MathML(node, all_species, inputs, level=0, logic_function=None):
    if logic_function is None:
        logic_function = []

    node_name = node.getName()
    #     print("  "*level, node_name)

    if node.getNumChildren() == 0:

        if node_name in all_species:
            return node_name
        if node_name in inputs:
            return inputs[node_name]["threshold"]
        if node.isInteger():
            return node.getInteger()
        print(dir(node))

        raise ValueError(
            f"Unspecified input '{node_name}' in transition function!")

    children = []
    for i in range(node.getNumChildren()):
        child = node.getChild(i)
        children.append(parse_MathML(child, all_species, inputs, level + 1))

    operator = node_name

    if operator in ["and", "times"]:
        operator = "&"  # Treat "times" as a logical "and"
    elif operator in ["or", "plus"]:
        operator = "|"  # Treat "plus" as a logical "or"
    elif operator == "xor":
        return handle_xor(children)
    elif operator in ["eq", "neq", "gt", "lt", "geq", "leq"]:
        return handle_comparison(operator, children)
    elif operator == "not":
        return handle_not(children)
    else:
        raise ValueError(f"Unsupported math symbol: {operator}!")

    if level > 0:
        return "(" + f" {operator} ".join(children) + ")"
    return f" {operator} ".join(children)


################################################################################
### OUTPUT
################################################################################


def booldog2sbmlqual(network, outfile, level=3, version=1, qual_version=1):
    '''Transforms booldog::Network to model in SBML-qual format '''

    # use to make nice species ids
    alphnumeric_pattern = re.compile(r'[\W_]+')

    sbmlns = SBMLNamespaces(level, version, "qual", qual_version)
    sbmldoc = SBMLDocument(sbmlns)

    # mark qual as required
    sbmldoc.setPackageRequired("qual", True)

    # create the Model
    model = sbmldoc.createModel()

    # create the Compartment
    compartment = model.createCompartment()
    compartment.setId("c")
    compartment.setConstant(True)

    # Get a QualModelPlugin object plugged in the model object.
    mplugin = model.getPlugin("qual")

    node_dict = {}
    for i, node in enumerate(sorted(network.nodes)):
        # id_ = f"s{i}"
        id_ = alphnumeric_pattern.sub('', node.lower())
        qs = mplugin.createQualitativeSpecies()
        qs.setId(id_)
        qs.setCompartment("c")
        if network.is_constant(node):
            qs.setConstant(True)
        else:
            qs.setConstant(False)
        qs.setName(node)
        node_dict[node] = id_

    for i, line in enumerate(network.to_bnet(header=False).split("\n")):
        if line == "":
            continue

        # print("--")
        # print(line)
        # print("--")
        target, rule = [s.strip() for s in line.split(",", 1)]

        if rule == "":
            continue

        node = target.strip()

        id_ = f"t{i}"
        node_id = node_dict[node]

        # create the Transition
        t = mplugin.createTransition()
        t.setId(id_)

        input_node_dict = {}
        for input_node in sorted(network.get_parents(node)):
            input_species_id = node_dict[input_node]
            input_node_id = f"theta_{id_}_{input_species_id}"
            input_node_dict[input_node] = input_node_id

            t_input = t.createInput()
            t_input.setId(input_node_id)
            t_input.setQualitativeSpecies(input_species_id)
            t_input.setTransitionEffect(libsbml.INPUT_TRANSITION_EFFECT_NONE)
            t_input.setThresholdLevel(1)

        t_output = t.createOutput()
        t_output.setId(f"{id_}_{node_id}")
        t_output.setQualitativeSpecies(node_id)
        t_output.setTransitionEffect(
            libsbml.OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL)

        ft = t.createFunctionTerm()

        ft.setResultLevel(1)
        try:
            math = rule2formula(rule, node_dict, input_node_dict)
            ft.setMath(math)
        except Exception as err:
            raise BoolDogSBMLException(
                f"Fatal error occurred in preparing transition for line: '{line}'"
            ) from err

        dt = t.createDefaultTerm()
        dt.setResultLevel(0)

    num_errors = sbmldoc.getNumErrors()
    if num_errors > 0:
        logger.warning("Generated SBML file has %i error(s):", num_errors)
        for i in range(num_errors):
            logger.warning("  %s", sbmldoc.getError(i).getMessage())

    libsbml.writeSBML(sbmldoc, outfile)
    logger.info('Wrote Network as a Boolean model in SBML-qual to %s', outfile)


def rule2formula(rule, node_dict, input_node_dict):
    ''' Rule (as in bnet) to a formula as supported by parseL3Formula

    &  --> &&
    |  --> ||

    A  --> (A >= theta_A)
    !A --> (A < theta_A)

    '''

    new_factors = []
    for factor in rule.split("|"):
        nodes = [s.strip() for s in factor.split("&")]
        new_nodes = []
        for node in nodes:
            if node.startswith("!"):
                node = node.lstrip("!")
                node_id = node_dict[node]
                new_node = f"({node_id} < {input_node_dict[node]})"
            else:
                node_id = node_dict[node]
                new_node = f"({node_id} >= {input_node_dict[node]})"
            new_nodes.append(new_node)
        new_factors.append(" && ".join(new_nodes))
    new_rule = " || ".join(new_factors)

    math = libsbml.parseL3Formula(new_rule)

    return math
