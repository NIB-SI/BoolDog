# Known bugs

Originally found while doing a full documentation pass over `booldog/`
(2026-07-14). All fixed on 2026-07-15 except the items in the last section,
which were investigated and are not actually bugs.

## Fixed

- `booldog/utils/misc.py`, `file_writable`: opened the target file in
  `'wb'` mode to check writability, which truncated (or created) it as a
  side effect *before* the real write. This was actively harmful in two
  call sites: `BooleanSimulationResult.export()` (an unimplemented stub)
  called it and did nothing else, silently destroying any pre-existing
  file at that path for no benefit; and `_to_writer` (used by
  `to_bnet`/`to_primes`/`to_sbmlqual`) called it before a writer function
  that can fail partway through building its output (e.g. `write_sbmlqual`
  building an SBML document), permanently losing the original file's
  content even though the export never actually completed. Fixed by: (1)
  removing the stub's call entirely (it isn't implemented yet, so nothing
  needs to write there), and (2) changing `file_writable` to open an
  already-existing path in `'r+b'` mode (read+write, no truncation)
  instead of `'wb'`, so a successful check no longer destroys existing
  content; the not-yet-existing-path case is unchanged.

- `booldog/boolean/modifications.py`: `modify_network()` read
  `modification.modification_type` (nonexistent attribute — `Modification`
  only sets `self.type`) and called a nonexistent `self.update_rule(...)`
  for the "update" case. Fixed to read `.type` and call `update_node`.

- `booldog/simulation_result/boolean_result.py`:
  `BooleanStateSpace.plot_state_space` raised `NameError` whenever
  `plot_nodes` was given, because it filtered the undefined name `nodes`
  instead of the `plot_nodes` parameter. Fixed.

- `booldog/network.py`: `BoolDogModel.__init__` never set `self.modelinfo`
  when `modelinfo=None`. Now defaults to a plain `BoolDogModelInfo()`
  instance, so code that dereferences `.modelinfo` unconditionally
  (`booldog/io/cytoscape.py`, `booldog/simulation_result/continuous_result.py`)
  no longer raises `AttributeError`.

- `booldog/simulation_result/continuous_result.py`:
  `ContinuousSimulationResult.export()` unconditionally read
  `self.ode_system.param_dict`, which only existed on `BooleCubeODE`. Fixed
  by giving `SquadODE` its own `param_dict` (`{"gamma": ..., "h": ...}`),
  matching `BooleCubeODE`'s.

- `booldog/boolean/modifications.py:300` — `f"{node}` typo (the `f` was
  inside the string literal). Fixed to a real f-string, and to reference
  the correct variable (`node_id`).

- `booldog/continuous/semi_quantitative.py` — same class of bug:
  `"        {node} -> released\n"` wasn't an f-string. Fixed.

- `booldog/io/interaction_networks.py`, `read_networkx`: `node_name_key`
  was accepted but ignored (hardcoded `"name"`). Fixed to use the
  parameter.

- `booldog/io/interaction_networks.py`, `read_sif`: default
  `activator_symbol`/`inhibitor_symbol` were ints (`1`/`-1`), but SIF
  values are always parsed as strings, so the defaults never matched
  anything. `read_sif` now has its own string defaults (`"1"`/`"-1"`),
  verified against the SIF test fixture with no explicit override needed.

- `booldog/io/sbml.py`, `MathMLParser._handle_comparison`: constant/constant
  comparisons returned Python's `str(bool)` (`"True"`/`"False"`) instead of
  bnet's `"1"`/`"0"`. Fixed.

- `booldog/io/sbml.py`, `TransitionParser.parse_io`: an unrecognised
  species / unsupported transition effect / set output level on one
  output used `break`, dropping all subsequent outputs for that
  transition. Changed to `continue`, so only the offending output is
  skipped.

- `booldog/io/bnet.py`: `write_bnet(from_primes=True)` always returned
  `None`, discarding the string `pyboolnet.file_exchange.primes2bnet`
  returns when `outfile=None`. Fixed to return it, matching the
  `from_primes=False` contract.

- `booldog/io/interaction_logic.py`: `interactions2rules`'s return type
  hint was `Dict[str, Callable]`; fixed to `Dict[str, str]` (it returns
  bnet-format rule strings).

- `booldog/__init__.py`: `assert sys.version_info >= (3, 10)` didn't match
  the project's actual `>=3.12` requirement. Fixed.

- `booldog/continuous/ode_factory.py`: removed the dead `'placeholder'`
  transform code path (`ODE.__init__` returned before setting
  `self.boolean_network`, so using it would have crashed one line later;
  it was only reachable from already-commented-out code, itself now also
  removed from `docs/source/conf.py`).

- `booldog/continuous/ode_factory.py`, `SquadODE._get_system`:
  `off_nodes=[]` mutable default argument. Changed to `off_nodes=None`
  with the list built inside the function, matching `BooleCubeODE`'s
  existing pattern.

- `booldog/io/__init__.py`: `_to_writer` was dead code — every `to_x`
  method called its writer directly. `to_bnet`/`to_primes`/`to_sbmlqual`
  now route through it (verified with both positional and keyword
  `outfile` for all three).

- `booldog/classes.py`: removed the commented-out, fully unused
  `BoolDogRule` class.

## Investigated, not bugs

- `booldog/io/biomodels.py`'s double-slash URL
  (`f"{BIOMODELS_BASE_URL}/{model_id}..."` where `BIOMODELS_BASE_URL`
  already ends in `/`): the recorded VCR cassettes in `tests/cassettes/`
  show this exact double-slash URL was sent to the *real* BioModels API
  and got a successful response — the server tolerates it. Left as-is;
  "fixing" it would only risk breaking cassette matching for no benefit.

- `booldog/io/cytoscape.py`, `silence_p4c_loggers`: the logger names
  `"py4..."` and `"py4...S"` looked like truncated placeholder text, but
  checking py4cytoscape's own source
  (`py4cytoscape/py4cytoscape_logger.py`) confirms these are its actual,
  if oddly-chosen, real logger names (`detail_logger`/`summary_logger`).
  Verified `logging.getLogger("py4...")` returns the exact same `Logger`
  instance py4cytoscape logs through. Not a bug.

- `booldog/boolean/boolean.py`, `get_interactions`'s `# TODO this is
  wrong?` comment: independently re-verified both `direction='out'` and
  `direction='in'` against a network with a known regulator relationship
  (`node_C = node_A & node_B`) — both directions produce exactly the
  documented semantics. Comment removed.

- `booldog/utils/logger.py`, `setup_logger`: the commented-out
  `ch.setLevel(level)` is redundant, not missing — the logger's own level
  is already set (`logger.setLevel(level)`), and with a single handler
  left at `NOTSET`, that alone is sufficient to control verbosity.
  Removed the dead comment; no behaviour change.

- `booldog/utils/misc.py`, `parameter_to_array`: silently falls back to an
  all-ones array (with a logged warning) for input that isn't an
  int/float/dict/matching-length ndarray, rather than raising. Left
  as-is — this reads as a deliberate permissive-with-warning design
  choice for a parameter-coercion helper, not a bug, and changing it to
  raise would be a behaviour change rather than a fix.

