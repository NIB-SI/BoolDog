'''Benchmark a single bnet model, read from stdin: model load and
instantiation (prime-implicant conversion), Boolean simulation (from a
fixed initial state, i.e. the reachable subgraph only, not the full 2^n
state space), and continuous simulation.

Used identically for both benchmark axes (real models and synthetic
in-degree-controlled models, see models.py) -- run as a subprocess (not
imported) so the driver script (run_benchmark.py) can bound each data
point with a wall-clock timeout independently, without a pathologically
slow high-in-degree model (see KNOWN_BUGS.md / hklarner/pyboolnet#122)
blocking the rest of the batch.
'''
import json
import sys
import time

import booldog


def benchmark(bnet_text, t_max):
    row = {}

    t0 = time.perf_counter()
    bn = booldog.BoolDogModel.from_bnet(bnet_text)
    row["n_nodes"] = len(bn.nodes)
    row["import_and_primes_s"] = time.perf_counter() - t0

    init = {nid: 0 for nid in bn.node_ids}
    t0 = time.perf_counter()
    bn.boolean_simulation(initial_states=init)
    row["boolean_sim_s"] = time.perf_counter() - t0

    t0 = time.perf_counter()
    bn.continuous_simulation(t_max=t_max, transform="normalisedhillcube")
    row["continuous_sim_s"] = time.perf_counter() - t0

    return row


if __name__ == "__main__":
    print(json.dumps(benchmark(sys.stdin.read(), float(sys.argv[1]))))
