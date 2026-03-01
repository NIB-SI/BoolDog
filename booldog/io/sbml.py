'''Read-write SBML-qual files

Notes
=====

SBML-qual code is inspired by:

* loadSBML.R of BoolNet (January 2024)
    https://github.com/cran/BoolNet/blob/af3a714c5bfa72ee7507db9c4eaf90ba2cd91809/R/loadSBML.R
* sbmlqual.py of CellNOpt (January 2024)
    https://github.com/cellnopt/cellnopt/blob/319c956e890ed7ae88c8e6e540de7061acfaae00/cno/io/sbmlqual.py
'''

# The function `parse_mathml` (and its dependencies) was initially written by
# ChatGPT, and then fixed by hand.

import re
import logging
from pathlib import Path

try:
    import libsbml
    from libsbml import SBMLReader, SBMLNamespaces, SBMLDocument
    _SBML_AVAILABLE = True
except ImportError as e:
    _SBML_AVAILABLE = False

from booldog.utils.boolean_normal_forms import functions2mindnf
from booldog.classes import BoolDogNode, BoolDogModelInfo

logger = logging.getLogger(__name__)

# libsbml sets certain missing values to SBML_MAX_INT, so to determine if they
# are missing, need to check against SBML_MAX_INT. Value from libsbml C code.
# TODO: is there a better option?
SBML_INT_MAX = 2147483647

# BoolNet rule tokens
TOKEN_REGEX = re.compile(r"""
    \s*
    (
        \!  |         # !
        \&  |         # &
        \|  |         # |
        \(  |         # (
        \)  |         # )
        [A-Za-z_][A-Za-z0-9_]*  # identifier (letters, numbers, underscore)
    )
""", re.VERBOSE)

##############################
# Utility classes/functions
##############################


class BoolDogSBMLException(Exception):
    '''Custom Exception for SBML parsing '''


class SBMLQualReader:

    def __init__(self, file):
        '''Read SBML-qual xml file to a BooleanNetwork object, via bnet.

        Parameters
        ----------
        file : str or path-like
            Path to SBML-qual file containing a Boolean network.

        '''

        if not _SBML_AVAILABLE:
            raise ImportError("libsbml is not available.")
        self.file = file
        self.document = SBMLReader().readSBML(file)

        if self.document.getNumErrors() > 0:
            for i in range(self.document.getNumErrors()):
                logger.error("Error: %s",
                             self.document.getError(i).getMessage())
            raise BoolDogSBMLException("SBML file contains errors.")

        self.model = self.document.getModel()
        self.model_id = self.model.getId()

        self.plugin = self._get_qual_plugin()

        species = self._get_all_species()
        self.all_species = [s[0] for s in species]
        self.species_names = {s[0]: s[1] for s in species}

        self.transitions = self._get_all_transitions()
        self.rules = self._get_all_rules()

    def _get_qual_plugin(self):
        for i in range(self.document.getNumPlugins()):
            plugin = self.model.getPlugin(i)
            if plugin.getPackageName() == 'qual':
                return plugin
        raise BoolDogSBMLException("SBML file missing 'qual' plugin.")

    def _get_all_species(self):
        return [(self.plugin.getQualitativeSpecies(i).getId(),
                 self.plugin.getQualitativeSpecies(i).getName())
                for i in range(self.plugin.getNumQualitativeSpecies())]

    def _get_all_transitions(self):
        return [
            self.plugin.getTransition(i)
            for i in range(self.plugin.getNumTransitions())
        ]

    def _get_all_rules(self):

        rules = {}
        for transition in self.transitions:
            inputs, outputs = TransitionParser.parse_io(
                transition, self.all_species)
            rule = TransitionParser.parse_function(transition,
                                                   self.all_species, inputs)
            for output in outputs:
                rules[output['species']] = rule

                logger.debug("Final rule: %s : %s", output['species'], rule)

        return rules

    def to_bnet(self):
        '''Converts the SBML-qual file to a bnet format.

        Returns
        -------
        bnet: str
            bnet representation of the Boolean network.
        '''

        bnet = "targets, factors\n"
        for target, rule in self.rules.items():
            bnet += f"{target}, {rule}\n"
        return bnet


class TransitionParser:
    '''Parse SBML transition to bnet format'''

    @staticmethod
    def parse_io(transition, all_species):
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
                "id":
                sp.getId(),
                "species":
                sp.getQualitativeSpecies(),
                "sign":
                sp.getSign(),
                "threshold":
                None if sp.getThresholdLevel() == SBML_INT_MAX else
                sp.getThresholdLevel(),
                "transition_effect":
                sp.getTransitionEffect()
            }

            if d["species"] not in all_species:
                logger.warning("Species '%s' not defined in model",
                               d["species"])

            if d["transition_effect"] != libsbml.INPUT_TRANSITION_EFFECT_NONE:
                logger.warning(
                    "Effect '%s' is not 'None' for transition input %s",
                    d["transition_effect"], d["id"])

            if not (d["threshold"] in (0, 1)):
                logger.warning(
                    "Threshold '%s' is not 0 or 1 for transition input %s",
                    d["threshold"], d["id"])

            inputs[d["id"]] = d

        outputs = []
        for sp in transition.getListOfOutputs():
            d = {
                "species":
                sp.getQualitativeSpecies(),
                "transition_effect":
                sp.getTransitionEffect(),
                "output_level":
                None if sp.getOutputLevel() == SBML_INT_MAX else
                sp.getOutputLevel()
            }

            if not (d["species"] in all_species):
                logger.warning("Species '%s' not defined in model",
                               d["species"])
                break

            if d["transition_effect"] != libsbml.OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL:
                logger.warning(
                    "Transition effect '%s' not defined for Boolean model",
                    d["transition_effect"])
                break

            if d["output_level"] is not None:
                logger.warning(
                    "Output level '%s' not supported for Boolean model",
                    d["output_level"])
                break
            outputs.append(d)

        return inputs, outputs

    @staticmethod
    def parse_function(transition, all_species, inputs):
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
        activation, inhibition = [], []

        for term in function_terms:

            output_level = term.getResultLevel()
            if output_level not in (0, 1):
                logger.warning('ResultLevel %i is not 0 or 1', output_level)

            logger.debug('ResultLevel %i', output_level)

            if term.getName().lower() == "defaultterm":
                continue

            try:
                rule = MathMLParser.parse(term.getMath(), all_species, inputs)
            except Exception as err:
                raise BoolDogSBMLException(
                    f"Failed parsing transition {transition.getId()}") from err

            logger.debug('%s <==> %s',
                         libsbml.formulaToL3String(term.getMath()), rule)

            if output_level == 1:
                activation.append(rule)
            else:
                inhibition.append(rule)

        default_term = transition.getDefaultTerm().getResultLevel()

        # TODO TODO TODO!!!
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
        logger.debug("Activation: %s", activation)
        logger.debug("Inhibition: %s", inhibition)
        if len(activation) > 0:
            rule = f"( {' | '.join(activation)} )"
            if len(inhibition) > 0:
                rule += f" & !( {' | '.join(inhibition)} )"
        elif len(inhibition) > 0:
            rule = f"!( {' | '.join(inhibition)} )"
        else:
            rule = str(default_term)

        return rule


class MathMLParser:
    '''Parse MathML to bnet format'''

    @staticmethod
    def parse(node, all_species, inputs, level=0):
        node_name = node.getName()
        if node.getNumChildren() == 0:
            if node_name in all_species:
                return node_name
            if node_name in inputs:
                return inputs[node_name]["threshold"]
            if node.isInteger():
                return node.getInteger()
            raise ValueError(
                f"Unspecified input '{node_name}' in transition function!")

        children = [
            MathMLParser.parse(node.getChild(i),
                               all_species,
                               inputs,
                               level=level + 1)
            for i in range(node.getNumChildren())
        ]

        return MathMLParser._handle_operator(node_name, children, level)

    @staticmethod
    def _handle_operator(operator, children, level):

        if operator in ["and", "times"]:
            op = " & "  # Treat "times" as a logical "and"
        elif operator in ["or", "plus"]:
            op = " | "  # Treat "plus" as a logical "or"
        elif operator == "xor":
            return MathMLParser._handle_xor(children)
        elif operator in ["eq", "neq", "gt", "lt", "geq", "leq"]:
            return MathMLParser._handle_comparison(operator, children)
        elif operator == "not":
            if len(children) != 1:
                raise ValueError(
                    "Unary operator \"not\" can only have one argument!")
            return f"!{children[0]}"
        else:
            raise ValueError(f"Unsupported math operator: {operator}")

        return f"({op.join(children)})" if level > 0 else f"{op.join(children)}"

    @staticmethod
    def _handle_xor(children):
        '''Create disjunctive normal form (DNF) of an xor '''

        def xor(*l):
            return not (sum(l) % 2 == 0)

        xor.depends = children
        return functions2mindnf({"xor_func": xor})["xor_func"]

    @staticmethod
    def _handle_comparison(operator, children):
        '''Return bnet form of an (mathml) operator between two children.'''

        if len(children) != 2:
            raise ValueError(f"{operator} requires two operands.")

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


class SBMLQualWriter:

    def __init__(self, network, level=3, version=1, qual_version=1):
        self.network = network

        ns = SBMLNamespaces(level, version, "qual", qual_version)
        doc = SBMLDocument(ns)

        # mark qual as required
        doc.setPackageRequired("qual", True)

        # create model
        model = doc.createModel()

        compartment = model.createCompartment()
        compartment.setId("c")
        compartment.setConstant(True)

        # get a QualModelPlugin object plugged in the model object.
        self.mplugin = model.getPlugin("qual")

        self.node_dict = {}

        self._add_species()
        self._add_transitions()

        if num_errors := doc.getNumErrors() > 0:
            logger.warning("Generated SBML file has %i error(s):", num_errors)
            for i in range(num_errors):
                logger.warning("SBML error: %s", doc.getError(i).getMessage())

        self.doc = doc

    def write(self, outfile):

        if isinstance(outfile, Path):
            outfile = str(outfile)

        libsbml.writeSBML(self.doc, outfile)
        logger.info('Wrote Network as a Boolean model in SBML-qual to %s',
                    outfile)

    def _add_species(self):
        for node_id, node in self.network.nodes.items():
            node_id_sbml = re.sub(r'[\W_]+', '', node_id.lower())
            species = self.mplugin.createQualitativeSpecies()
            species.setId(node_id_sbml)
            species.setCompartment("c")
            species.setConstant(self.network.is_constant(node))
            species.setName(node.name)
            self.node_dict[node_id] = node_id_sbml

    def _add_transitions(self):
        for target in self.network.nodes.values():
            rule = target.rule
            target_id = self.node_dict[target.identifier]

            transition_id = f"tr_{target_id}"
            transition = self.mplugin.createTransition()
            transition.setId(transition_id)

            # Needed for future implementation of thresholds
            input_node_dict = {} # local to this transition, links input node to threshold
            # for parent in sorted(self.network.get_parents(target)):
            #     parent_id = node_dict[parent]
            #     theta_id = f"{transition_id}_{parent_id}_theta"
            #     input_node_dict[parent] = theta_id
            #     inp = transition.createInput()
            #     inp.setId(theta_id)
            #     inp.setQualitativeSpecies(parent_id)
            #     inp.setTransitionEffect(libsbml.INPUT_TRANSITION_EFFECT_NONE)
            #     inp.setThresholdLevel(1)

            out = transition.createOutput()
            out.setId(f"{transition_id}_out")
            out.setQualitativeSpecies(target_id)
            out.setTransitionEffect(
                libsbml.OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL)

            func = transition.createFunctionTerm()
            func.setResultLevel(1)

            try:
                func.setMath(
                    self._rule_to_formula(rule, input_node_dict))
            except Exception as err:
                raise BoolDogSBMLException(
                    f"Fatal error occurred in preparing transition for {target}: '{rule}'"
                ) from err

            transition.createDefaultTerm().setResultLevel(0)

    def _rule_to_formula(self, rule, input_node_dict):
        ''' Rule (as in bnet) to a formula as supported by parseL3Formula

        Notes
        -----

        Replaces:

        * `&`  with `&&`
        * `|`  with `||`

        To support thresholds, we need to replace the bnet format with a
        format supported by libsbml's parseL3Formula:

        * `A`  to `(A >= theta_A)` f"({node_dict[node]} >= {input_node_dict[node]})"
        * `!A` to `(A < theta_A)`  f"({node_dict[node]} < {input_node_dict[node]})"

        This is currently NOT implemented.
        '''


        tokens = TOKEN_REGEX.findall(rule)
        output = []

        for tok in tokens:

            # --- logical operators ---
            if tok == "&":
                output.append(" && ")
            elif tok == "|":
                output.append(" || ")

            # preserve existing tokens
            elif tok in ("!", "(", ")"):
                output.append(tok)

            # --- identifier ---
            else:
                output.append(self.node_dict[tok])

        return libsbml.parseL3Formula("".join(output))

###############################
# In
###############################

def read_sbmlqual(file):
    ''' Create Network from a SBML-qual file

    Parameters
    ----------
    file : str
        Path to SBML-qual file.

    Notes
    -----

    The SBML-qual file is converted to a Boolean network using libsbml, via
    the bnet format. To access the bnet format directly, you can use `py:booldog.io.sbmlqual2bnet`.


    '''
    if not _SBML_AVAILABLE:
        raise ImportError(
            'libsbml (https://sbml.org/software/libsbml/libsbml-docs/api/python/) '
            'is needed to read models in SBML format. '
            'We suggest you install it using pip. ')

    reader = SBMLQualReader(file)
    nodes = []
    for node_id in reader.all_species:
        if node_id in reader.rules:
            rule = reader.rules[node_id]
        else:
            rule = None
        node = BoolDogNode(
            identifier=node_id,
            rule=rule,

            # additional info from SBML, for now just the name,
            # but could be extended to include other attributes
            # (e.g., initial value, compartment, annotations, etc.)
            name=reader.species_names[node_id])
        nodes.append(node)

    modelinfo = BoolDogModelInfo(identifier=reader.model_id,
                                 source=file,
                                 source_format="sbml-qual")

    return {"nodes": nodes, "modelinfo": modelinfo, "primes": None}


################################
# Out
################################


def write_sbmlqual(model, outfile, **kwargs):
    ''' '''

    # IDEA: Keep the original (and updated) logic rules as
    # attribute, and have an option to use them to create the sbml -qual
    # transitions, instead of using the bnet/primes DNF

    if _SBML_AVAILABLE:
        writer = SBMLQualWriter(model)
        writer.write(outfile, **kwargs)
        return

    raise ImportError(
        'libsbml (https://sbml.org/software/libsbml/libsbml-docs/api/python/) '
        'is needed to write models to SBML format. '
        'We suggest you install it using pip. ')
