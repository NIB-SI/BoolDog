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
'''int: Sentinel value libsbml uses for "unset" integer attributes (e.g. a
qualitative species/transition's threshold or output level). Compared
against to determine whether such an attribute was actually set in the
SBML file.'''

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
'''re.Pattern: Tokenizes a bnet-format rule string into ``!``, ``&``, ``|``,
``(``, ``)`` and identifier tokens, ignoring surrounding whitespace. Used by
:meth:`SBMLQualWriter._rule_to_formula` to translate bnet rule syntax into
the ``&&``/``||`` syntax accepted by ``libsbml.parseL3Formula``.'''

##############################
# Utility classes/functions
##############################


class BoolDogSBMLException(Exception):
    '''Custom Exception for SBML parsing '''


class SBMLQualReader:
    '''Read an SBML-qual XML file, converting each transition to a
    bnet-format Boolean rule.

    Parameters
    ----------
    file : str or path-like
        Path to SBML-qual file containing a Boolean network.

    Raises
    ------
    ImportError
        If libsbml is not installed.
    BoolDogSBMLException
        If the SBML document fails to parse without errors, has no 'qual'
        plugin, or a transition's function term fails to parse (see
        :meth:`TransitionParser.parse_function`).
    '''

    def __init__(self, file):

        if not _SBML_AVAILABLE:
            raise ImportError("libsbml is not available.")
        self.file = file
        '''str or path-like: Path to the SBML-qual file, as given at
        construction.'''
        self.document = SBMLReader().readSBML(file)
        '''libsbml.SBMLDocument: The parsed SBML document.'''

        if self.document.getNumErrors() > 0:
            for i in range(self.document.getNumErrors()):
                logger.error("Error: %s",
                             self.document.getError(i).getMessage())
            raise BoolDogSBMLException("SBML file contains errors.")

        self.model = self.document.getModel()
        '''libsbml.Model: The SBML model contained in :attr:`document`.'''
        self.model_id = self.model.getId()
        '''str: The SBML model's id attribute.'''

        self.plugin = self._get_qual_plugin()
        '''libsbml.QualModelPlugin: The 'qual' package plugin of
        :attr:`model`, providing access to qualitative species and
        transitions.'''

        species = self._get_all_species()
        self.all_species = [s[0] for s in species]
        '''list of str: Ids of all qualitative species defined in the
        model.'''
        self.species_names = {s[0]: s[1] for s in species}
        '''dict: Mapping of qualitative species id to its (possibly empty)
        SBML name attribute.'''

        self.transitions = self._get_all_transitions()
        '''list of libsbml.Transition: All transitions defined in the
        model.'''
        self.rules = self._get_all_rules()
        '''dict: Mapping of output species id to its bnet-format Boolean
        rule string, derived from the model's transitions (see
        :meth:`_get_all_rules`).'''

    def _get_qual_plugin(self):
        '''Find and return the model's 'qual' package plugin.

        Returns
        -------
        libsbml.QualModelPlugin
            The 'qual' plugin attached to :attr:`model`.

        Raises
        ------
        BoolDogSBMLException
            If the model has no plugin with package name 'qual'.
        '''
        for i in range(self.document.getNumPlugins()):
            plugin = self.model.getPlugin(i)
            if plugin.getPackageName() == 'qual':
                return plugin
        raise BoolDogSBMLException("SBML file missing 'qual' plugin.")

    def _get_all_species(self):
        '''List all qualitative species defined via the 'qual' plugin.

        Returns
        -------
        list of tuple
            ``(id, name)`` pairs for every qualitative species in
            :attr:`plugin`, in declaration order.
        '''
        return [(self.plugin.getQualitativeSpecies(i).getId(),
                 self.plugin.getQualitativeSpecies(i).getName())
                for i in range(self.plugin.getNumQualitativeSpecies())]

    def _get_all_transitions(self):
        '''List all transitions defined via the 'qual' plugin.

        Returns
        -------
        list of libsbml.Transition
            Every transition in :attr:`plugin`, in declaration order.
        '''
        return [
            self.plugin.getTransition(i)
            for i in range(self.plugin.getNumTransitions())
        ]

    def _get_all_rules(self):
        '''Derive a bnet-format Boolean rule for each transition's output
        species.

        For each transition, inputs/outputs are extracted with
        :meth:`TransitionParser.parse_io` and its function term is parsed
        with :meth:`TransitionParser.parse_function`; the resulting rule
        string is assigned to every one of that transition's output
        species.

        Returns
        -------
        dict
            Mapping of output species id to its bnet-format rule string.
        '''

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
            Dictionary mapping each input's id to a dict of its species
            information: ``"id"``, ``"species"`` (qualitative species id),
            ``"sign"``, ``"threshold"`` (int, or None if unset in the SBML
            file), and ``"transition_effect"``.
        outputs: list
            List of dicts (one per transition output), each with keys
            ``"species"``, ``"transition_effect"`` and ``"output_level"``
            (int, or None if unset).

        Notes
        -----
        Warnings are logged (not raised) for inputs/outputs referencing a
        species not in `all_species`, inputs with a transition effect other
        than "None", or inputs with a threshold that isn't 0 or 1 — these
        are conditions this Boolean-only SBML-qual reader does not support,
        but does not treat as fatal for inputs.

        For outputs, an unrecognised species, an unsupported transition
        effect (i.e. not "assignmentLevel"), or a set output level are all
        logged as warnings and skipped (excluded from the returned
        `outputs` list), without affecting collection of any other
        outputs for this transition.
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
                continue

            if d["transition_effect"] != libsbml.OUTPUT_TRANSITION_EFFECT_ASSIGNMENT_LEVEL:
                logger.warning(
                    "Transition effect '%s' not defined for Boolean model",
                    d["transition_effect"])
                continue

            if d["output_level"] is not None:
                logger.warning(
                    "Output level '%s' not supported for Boolean model",
                    d["output_level"])
                continue
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

        Raises
        ------
        BoolDogSBMLException
            If :meth:`MathMLParser.parse` fails to parse a function term's
            MathML.

        Notes
        -----
        Function terms named "defaultTerm" are skipped (the actual default
        result is read separately from ``transition.getDefaultTerm()``).
        Every other function term is parsed to a bnet-format expression via
        :meth:`MathMLParser.parse`, and sorted by its ``ResultLevel`` into
        an ``activation`` list (level 1) or ``inhibition`` list (level 0 or
        anything else).

        The final rule is built as:

        * ``"( a1 | a2 | ... ) & !( i1 | i2 | ... )"`` if there is at least
          one activation and one inhibition term;
        * ``"( a1 | a2 | ... )"`` if there are only activation terms;
        * ``"!( i1 | i2 | ... )"`` if there are only inhibition terms;
        * ``str(default_term)`` (i.e. ``"0"`` or ``"1"``) if there are no
          non-default function terms at all.
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
    '''Recursively parse a libsbml MathML AST (as used in SBML-qual
    FunctionTerms) into a bnet-format Boolean rule string.'''

    @staticmethod
    def parse(node, all_species, inputs, level=0):
        '''Recursively parse a MathML AST node to bnet syntax.

        Parameters
        ----------
        node : libsbml.ASTNode
            The (sub-)expression to parse.
        all_species : list
            List of ids of all species present in the model; leaf nodes
            named after one of these are treated as species references.
        inputs : dict
            Mapping of transition input id to its species information (as
            returned by :meth:`TransitionParser.parse_io`); leaf nodes
            named after one of these keys are resolved to that input's
            ``"threshold"`` value.
        level : int, optional
            Current recursion depth, used by :meth:`_handle_operator` to
            decide whether to wrap the result in parentheses (top-level
            expressions, ``level == 0``, are not parenthesised). Default 0.

        Returns
        -------
        str or int
            For a leaf node: the species id (str) if it names a species in
            `all_species`; the referenced input's threshold (int or None)
            if it names a key in `inputs`; or the node's integer value
            (int) if it is an integer literal. For an internal node: the
            bnet-format string built by :meth:`_handle_operator` from its
            (recursively parsed) children.

        Raises
        ------
        ValueError
            If a leaf node's name matches neither a species nor an input
            id, and it is not an integer literal.
        '''
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
        '''Combine already-parsed child expressions with a MathML
        operator, into a bnet-format expression string.

        Parameters
        ----------
        operator : str
            MathML operator name (e.g. "and", "or", "times", "plus", "xor",
            "not", or a comparison: "eq", "neq", "gt", "lt", "geq", "leq").
        children : list
            The operator's operands, already parsed to bnet-format strings
            (or, for comparisons, possibly ints — see
            :meth:`_handle_comparison`).
        level : int
            Recursion depth of the containing expression (see
            :meth:`parse`); if greater than 0 the combined "and"/"or"
            expression is wrapped in parentheses.

        Returns
        -------
        str
            The combined bnet-format expression. "times"/"plus" are
            treated as "and"/"or" respectively (a common encoding when
            Boolean values are represented as 0/1 numerically).

        Raises
        ------
        ValueError
            If `operator` is "not" and does not have exactly one child, or
            if `operator` is not one of the recognised MathML operators.
        '''

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
        '''Combine already-parsed child expressions with an n-ary xor,
        via a minimal disjunctive normal form (DNF).

        Since bnet syntax has no native xor operator, each of `children`
        (bnet-format expression strings) is treated as an opaque Boolean
        variable, an n-input xor truth table over them is built, and
        :func:`booldog.utils.boolean_normal_forms.functions2mindnf` is used
        to minimise it to a DNF bnet-format expression string, substituting
        the original child expression strings back in as "variable names".

        Parameters
        ----------
        children : list of str
            The xor's operands, already parsed to bnet-format strings.

        Returns
        -------
        str
            A minimal-DNF bnet-format expression equivalent to the xor of
            `children`.
        '''

        def xor(*l):
            return not (sum(l) % 2 == 0)

        xor.depends = children
        return functions2mindnf({"xor_func": xor})["xor_func"]

    @staticmethod
    def _handle_comparison(operator, children):
        '''Return the bnet form of a MathML comparison operator applied to
        two already-parsed operands.

        Parameters
        ----------
        operator : str
            One of "eq", "neq", "gt", "lt", "geq", "leq".
        children : list
            Exactly two operands: each either a bnet-format expression
            string (a Boolean-valued sub-expression) or an int (a
            threshold/integer-literal leaf value, from :meth:`parse`).

        Returns
        -------
        str
            The bnet-format result of applying `operator` to the two
            operands:

            * if both operands are int: the Python comparison's Boolean
              result, as bnet's ``"1"``/``"0"``;
            * if exactly one operand is an int (0 or 1): a simplified
              bnet-format expression referencing only the variable
              operand (e.g. ``x >= 1`` is just ``x``; ``x >= 0`` is
              always true, i.e. ``"1"``);
            * if neither operand is an int: a full bnet-format Boolean
              expression combining both operands' expression strings.

        Raises
        ------
        ValueError
            If `children` does not have exactly two elements.
        '''

        if len(children) != 2:
            raise ValueError(f"{operator} requires two operands.")

        is_const = [isinstance(child, int) for child in children]

        if all(is_const):
            children = [int(child) for child in children]
            bool_result = {
                "eq": lambda x, y: x == y,
                "neq": lambda x, y: x != y,
                "gt": lambda x, y: x > y,
                "lt": lambda x, y: x < y,
                "geq": lambda x, y: x >= y,
                "leq": lambda x, y: x <= y,
            }[operator](children[0], children[1])
            return "1" if bool_result else "0"

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
    '''Build and write an SBML-qual representation of a Boolean
    :class:`BoolDogModel`.

    Each node becomes a qualitative species (id sanitised to alphanumerics
    only), and each node's rule becomes a Transition with a single
    FunctionTerm (result level 1, built from the rule via
    :meth:`_rule_to_formula`) and a DefaultTerm of result level 0.

    Parameters
    ----------
    network : BoolDogModel
        The Boolean network to export.
    level : int, optional
        SBML level. Default 3.
    version : int, optional
        SBML version. Default 1.
    qual_version : int, optional
        Version of the 'qual' package. Default 1.
    '''

    def __init__(self, network, level=3, version=1, qual_version=1):
        self.network = network
        '''BoolDogModel: The network being exported, as given at
        construction.'''

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
        '''libsbml.QualModelPlugin: The 'qual' plugin of the SBML model
        being built, used to create qualitative species and transitions.'''

        self.node_dict = {}
        '''dict: Mapping of BoolDog node identifier to its sanitised
        SBML-qual species id, populated by :meth:`_add_species`.'''

        self._add_species()
        self._add_transitions()

        if num_errors := doc.getNumErrors() > 0:
            logger.warning("Generated SBML file has %i error(s):", num_errors)
            for i in range(num_errors):
                logger.warning("SBML error: %s", doc.getError(i).getMessage())

        self.doc = doc
        '''libsbml.SBMLDocument: The fully-built SBML document, ready to be
        written out with :meth:`write`.'''

    def write(self, outfile):
        '''Write the built SBML document to file.

        Parameters
        ----------
        outfile : str or Path
            Path to write the SBML-qual XML file to.

        Returns
        -------
        None
        '''

        if isinstance(outfile, Path):
            outfile = str(outfile)

        libsbml.writeSBML(self.doc, outfile)
        logger.info('Wrote Network as a Boolean model in SBML-qual to %s',
                    outfile)

    def _add_species(self):
        '''Create a qualitative species for every node in :attr:`network`.

        The SBML id is derived from the node identifier by lower-casing it
        and stripping any non-alphanumeric/underscore characters (so it is
        recorded in :attr:`node_dict`, since it may no longer match the
        original node identifier). Each species is placed in the single
        compartment "c", marked constant according to
        ``network.is_constant(node)``, and named after ``node.name``.

        Returns
        -------
        None
        '''
        for node_id, node in self.network.nodes.items():
            node_id_sbml = re.sub(r'[\W_]+', '', node_id.lower())
            species = self.mplugin.createQualitativeSpecies()
            species.setId(node_id_sbml)
            species.setCompartment("c")
            species.setConstant(self.network.is_constant(node))
            species.setName(node.name)
            self.node_dict[node_id] = node_id_sbml

    def _add_transitions(self):
        '''Create a Transition for every node in :attr:`network`, encoding
        its Boolean rule.

        Each transition has a single output (the node itself, with
        transition effect "assignmentLevel"), one FunctionTerm (result
        level 1) whose math is the node's rule converted via
        :meth:`_rule_to_formula`, and a DefaultTerm of result level 0.
        Transition inputs (regulator thresholds) are not currently created;
        see the commented-out code and :meth:`_rule_to_formula`'s Notes for
        the planned (not yet implemented) threshold support.

        Returns
        -------
        None

        Raises
        ------
        BoolDogSBMLException
            If building the MathML formula for a node's rule fails (wraps
            the underlying exception).
        '''
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
        ''' Convert a bnet-format rule to a MathML AST, as supported by
        ``libsbml.parseL3Formula``.

        Parameters
        ----------
        rule : str
            A Boolean rule in bnet format (using node identifiers as they
            appear in :attr:`network`, e.g. ``"A & !B"``).
        input_node_dict : dict
            Currently unused by this method (reserved/accepted for the
            planned threshold support described below, but not read here).

        Returns
        -------
        libsbml.ASTNode
            The MathML AST parsed from the reformatted formula (the return
            value of ``libsbml.parseL3Formula``), suitable for
            ``FunctionTerm.setMath``.

        Raises
        ------
        KeyError
            If an identifier token in `rule` is not a key of
            :attr:`node_dict` (i.e. not a node of :attr:`network`).

        Notes
        -----

        Tokenizes `rule` with :data:`TOKEN_REGEX` and, for each token,
        replaces:

        * `&`  with `&&`
        * `|`  with `||`
        * identifiers with their sanitised SBML species id
          (``self.node_dict[tok]``)
        * `!`, `(`, `)` are kept as-is

        then parses the reassembled string with
        ``libsbml.parseL3Formula``.

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
    ''' Parse an SBML-qual file into the data needed to construct a
    :py:class:`BoolDogModel`.

    Parameters
    ----------
    file : str
        Path to SBML-qual file.

    Returns
    -------
    data : dict
        Dictionary with keys ``"nodes"`` (list of :class:`BoolDogNode`, one
        per qualitative species; nodes with no associated transition get
        ``rule=None``), ``"modelinfo"`` (:class:`BoolDogModelInfo`, with
        ``identifier`` set to the SBML model id and ``source_format`` set
        to ``"sbml-qual"``), and ``"primes"`` (``None``, since primes are
        not computed by this reader). Suitable for ``BoolDogModel(**data)``.

    Raises
    ------
    ImportError
        If libsbml is not installed.

    Notes
    -----

    The SBML-qual file is converted to a Boolean network using libsbml, via
    the bnet format. To access the bnet format directly, you can construct
    a :class:`SBMLQualReader` and call its :meth:`SBMLQualReader.to_bnet`
    method.


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
    ''' Write a BoolDogModel object to an SBML-qual file, via
    :class:`SBMLQualWriter`.

    Parameters
    ----------
    model : BoolDogModel
        A BoolDog object representing a Boolean network.
    outfile : str or Path
        Path to write the SBML-qual XML file to.
    **kwargs
        Forwarded to :meth:`SBMLQualWriter.write`. ``SBMLQualWriter`` is
        also always constructed with its default
        ``level``/``version``/``qual_version``; those are not currently
        configurable from this function.

    Returns
    -------
    None

    Raises
    ------
    ImportError
        If libsbml is not installed.
    TypeError
        **Known bug, not intentional behaviour**: :meth:`SBMLQualWriter.write`
        only accepts ``outfile``, so passing any keyword arguments in
        ``**kwargs`` here raises `TypeError`. See ``KNOWN_BUGS.md``.
    '''

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
