# Known bugs

Found while doing a full documentation pass over `booldog/` (2026-07-14). Not
fixed — flagged here for follow-up since they were out of scope for a
docs-only change.

## Serious — likely broken features

- **`booldog/boolean/modifications.py`: `modify_network()` appears
  completely non-functional.** `Modification.__init__` only sets
  `self.type` (not `self.modification_type`), but `modify_network` reads
  `modification.modification_type` — this raises `AttributeError` before
  reaching any `case` branch. Separately, its `UPDATE` branch calls
  `self.update_rule(...)`, which doesn't exist anywhere in the codebase
  (the real method is `update_node`).

- **`booldog/simulation_result/boolean_result.py`:
  `BooleanStateSpace.plot_state_space` raises `NameError`.** When
  `plot_nodes` is truthy, it does
  `nodes = [node_id for node_id in nodes if node_id in self.network.node_ids]`
  — referencing `nodes` before it's ever assigned.

- **`booldog/network.py`: `BoolDogModel.__init__` never sets
  `self.modelinfo`** if `modelinfo=None` is passed (no fallback default).
  Downstream code dereferences `.modelinfo` unconditionally with no guard:
  `booldog/io/cytoscape.py` (`model.modelinfo.identifier`) and
  `booldog/simulation_result/continuous_result.py`
  (`self.network.modelinfo.source`) — both would raise `AttributeError`
  for a model constructed without `modelinfo=`.

- **`booldog/simulation_result/continuous_result.py`:
  `ContinuousSimulationResult.export()` unconditionally accesses
  `self.ode_system.param_dict`**, which only exists on `BooleCubeODE`, not
  `SquadODE` — breaks exporting any SQUAD-based simulation result.

## Smaller bugs

- `booldog/boolean/modifications.py:300` —
  `raise ValueError("f{node} is already present in Network.")`: the `f`
  prefix is typed *inside* the string literal, so it's never an f-string;
  the message literally reads `f{node} is already present...` instead of
  interpolating the node id.

- `booldog/continuous/semi_quantitative.py:~243` — same class of bug:
  `s += "        {node} -> released\n"` is a plain string, not an
  f-string, so `{node}` is never interpolated in that log message (the
  "starting a perturbation" branch just below it correctly uses an
  f-string).

- `booldog/io/interaction_networks.py`, `read_networkx`: the
  `node_name_key` parameter is accepted but never used — the code
  hardcodes `data.get("name", None)` instead of
  `data.get(node_name_key, None)`.

- `booldog/io/interaction_networks.py`, `read_graphml`: the
  `use_labels=True` parameter is accepted but never referenced anywhere in
  the function body — appears dead/unused.

- `booldog/io/interaction_networks.py`, `read_sif`: default
  `activator_symbol=1`/`inhibitor_symbol=-1` are ints, but SIF file values
  are always parsed as strings — the defaults likely never match unless
  the caller explicitly passes string symbols.

- `booldog/io/sbml.py`: `SBMLQualWriter.write(self, outfile)` doesn't
  accept `**kwargs`, but `write_sbmlqual` calls
  `writer.write(outfile, **kwargs)` — passing any kwargs raises
  `TypeError`.

- `booldog/io/sbml.py`, `MathMLParser._handle_comparison`: when both
  operands are constant ints, the result is `str(bool)`
  (`"True"`/`"False"`) rather than bnet's `"1"`/`"0"` — if joined into an
  `and`/`or` expression later, this would produce invalid bnet syntax.

- `booldog/io/sbml.py`, `TransitionParser.parse_io`: for a transition's
  outputs, an unrecognised species / unsupported transition effect / set
  output level all `break` out of the whole outputs loop rather than
  `continue`-ing past just that one output — a single malformed output
  silently drops any subsequent valid outputs for that transition.

- `booldog/io/bnet.py`: `write_bnet(model, outfile=None, from_primes=True)`
  always returns `None`, discarding the string
  `pyboolnet.file_exchange.primes2bnet` would otherwise return — this
  contradicts the "returns bnet string if outfile is None" contract that
  holds for `from_primes=False`.

- `booldog/io/biomodels.py`, `fetch_model`: the `sbml_file` parameter is
  unconditionally overwritten and has no effect.

- `booldog/io/biomodels.py`: the "first matching file" claim in the old
  docstring was wrong — the loop has no `break`, so it's actually the
  *last* matching file that's used.

- `booldog/io/biomodels.py`: `BIOMODELS_BASE_URL = "https://www.biomodels.org/"`
  combined with `f"{BIOMODELS_BASE_URL}/{model_id}..."` produces a double
  slash in the constructed request URL. Not verified against the live API
  — may or may not actually break requests.

- `booldog/io/interaction_logic.py`: `interactions2rules`'s return type
  hint is `Dict[str, Callable]`, but it actually returns
  `Dict[str, str]` (bnet-format rule strings, not callables).

- `booldog/io/circuit.py`, `booldog2circuit`: docstring claimed it returns
  `BooleanDiGraph`, but it actually returns a plain `nx.DiGraph(graph)`.

- `booldog/io/igraph.py` / `booldog/io/networkx.py`: both had a stale
  `as_logic_circuit` default in the docstring ("Default is False"); the
  actual signature default is `True`.

## Lower priority / worth a look

- `booldog/__init__.py:14` — `assert sys.version_info >= (3, 10)`, but
  `CLAUDE.md`/`pyproject.toml` both require Python `>=3.12`. Stale check.

- `booldog/continuous/ode_factory.py`: the `'placeholder'` transform path
  (`ode_factory(model, 'placeholder')`) looks dead/broken —
  `ODE.__init__` returns immediately without setting `self.boolean_network`
  when `transform == 'placeholder'`, so `BooleCubeODE.__init__` would
  crash on its next line. Only referenced from commented-out code.

- `booldog/continuous/ode_factory.py`, `SquadODE._get_system`:
  `off_nodes=[]` is a mutable default argument (classic Python gotcha) —
  currently only read, never mutated in place, so benign today.

- `booldog/io/__init__.py`: `_to_writer` is defined as the writer-side
  counterpart to `_from_reader`, but appears to be dead code — every
  `to_x` method calls its writer function directly instead of going
  through `_to_writer`.

- `booldog/io/cytoscape.py`, `silence_p4c_loggers`: the logger name
  strings passed look like truncated/placeholder text (e.g. `"py4..."`)
  rather than real py4cytoscape logger names — likely non-functional as
  written.

- `booldog/boolean/boolean.py`, `get_interactions`: has its own
  `# TODO this is wrong?` comment. Traced through `primes2igraph`'s edge
  semantics and the `direction='out'`/`'in'` branches appear to correctly
  implement the documented behavior, but this wasn't confirmed against
  any test coverage — worth a second look before removing the TODO.

- `booldog/utils/misc.py`, `parameter_to_array`: any `parameter` that
  isn't an int/float/dict/matching-length ndarray (including a
  wrong-length ndarray) silently falls through to a logged warning and
  returns an all-ones array, discarding the caller's value with no
  exception raised.

- `booldog/utils/misc.py`, `file_writable`: the writability "check" opens
  the target file in `'wb'` mode, which creates it (or truncates existing
  contents) as a side effect merely to test writability.

- `booldog/utils/logger.py`, `setup_logger`: `# ch.setLevel(level)` is
  commented out — unclear whether leaving the handler at `NOTSET`
  (deferring entirely to the logger's own level) is intentional or
  leftover debugging.

- `booldog/classes.py`: a commented-out `BoolDogRule` class — confirmed
  (via repo-wide grep) to be completely unused/dead code.
