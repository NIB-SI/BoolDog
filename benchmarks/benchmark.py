'''BoolDog performance/scalability benchmark: run the sweeps and/or plot
the results. See README.md for background and usage.

    python benchmark.py run    # write results/*.tsv
    python benchmark.py plot   # write results/figures/scaling.png
    python benchmark.py all    # both, in order

Two benchmark axes, both measuring the same three metrics (model load and
instantiation, Boolean simulation, continuous simulation; see
bench_worker.py) via the same generic sweep, each data point run as its
own subprocess with a wall-clock timeout so one pathologically slow
high-in-degree model doesn't block the rest of the batch (timeouts are
recorded, not silently dropped):

* node-count scaling, over real biological models (models.MODELS), with
  each model's actual max in-degree recorded alongside so runtime can be
  related to whichever variable actually drives it.
* in-degree scaling, over controlled synthetic single-node functions,
  isolating the mechanism (models.synthetic_indegree_bnet / KNOWN_BUGS.md).
'''
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from models import DATA_DIR, MODELS, fetch_all, max_in_degree, synthetic_indegree_bnet

HERE = Path(__file__).parent
RESULTS_DIR = HERE / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

TIMEOUT_S = 180
CONTINUOUS_T_MAX = 10
INDEGREE_VALUES = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

FIELDNAMES = ["n_nodes", "import_and_primes_s", "boolean_sim_s", "continuous_sim_s", "error"]

METRIC_LABELS = {
    "import_and_primes_s": "Model load and instantiation",
    "boolean_sim_s": "Boolean simulation (fixed initial state)",
    "continuous_sim_s": "Continuous simulation (t_max=10)",
}
MARKERS_BY_METRIC = {
    "Model load and instantiation": "o",
    "Boolean simulation (fixed initial state)": "X",
    "Continuous simulation (t_max=10)": "s",
}


# --------------------------------------------------------------------
# run
# --------------------------------------------------------------------

def run_sweep(items, stop_on_timeout=False):
    '''items: list of (extra_fields dict, bnet_text). Runs bench_worker.py
    on each bnet_text with a wall-clock timeout, returning one result row
    (extra_fields merged with the worker's JSON output, or an "error" key)
    per item.

    stop_on_timeout=True stops after the first timeout, on the assumption
    later items are at least as slow -- only valid when items are ordered
    by the actual variable driving runtime (e.g. the in-degree sweep).
    The node-count sweep deliberately mixes low- and high-in-degree models
    at similar sizes, so a timeout there says nothing about later items.
    '''
    rows = []
    for extra_fields, bnet_text in items:
        print(f"[{extra_fields}]...", file=sys.stderr)
        try:
            proc = subprocess.run(
                [sys.executable, str(HERE / "bench_worker.py"), str(CONTINUOUS_T_MAX)],
                input=bnet_text, capture_output=True, text=True, timeout=TIMEOUT_S,
            )
            # pyboolnet sometimes prints its own warnings directly to
            # stdout (not through logging), e.g. "WARNING The state
            # transition graph will consist of up to 2**24=... states" --
            # so pick out the one line that's actually our JSON, rather
            # than parsing all of stdout as JSON.
            result = None
            if proc.returncode == 0:
                for line in proc.stdout.splitlines():
                    if line.startswith("{"):
                        result = json.loads(line)
                        break
            error = None if result else (proc.stderr.strip().splitlines() or ["unknown error"])[-1]
        except subprocess.TimeoutExpired:
            result, error = None, "timeout"

        row = {**extra_fields, **(result or {"error": error})}
        print(f"  -> {row}", file=sys.stderr)
        rows.append(row)
        if error == "timeout" and stop_on_timeout:
            print(f"  (stopping sweep: timed out after {TIMEOUT_S}s)", file=sys.stderr)
            break
    return rows


def write_tsv(rows, fieldnames, path):
    RESULTS_DIR.mkdir(exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"wrote {path}")


def run_nodecount_benchmark():
    commit = fetch_all()
    items = []
    for model in MODELS:
        bnet_text = (DATA_DIR / f"{model['label']}.bnet").read_text()
        extra = {"label": model["label"], "dirname": model["dirname"],
                  "max_in_degree": max_in_degree(bnet_text)}
        items.append((extra, bnet_text))

    rows = run_sweep(items)
    write_tsv(rows, ["label", "dirname", "max_in_degree", *FIELDNAMES],
              RESULTS_DIR / "nodecount.tsv")

    (RESULTS_DIR / "biodivine_commit.txt").write_text(commit + "\n")
    booldog_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=HERE, capture_output=True, text=True, check=True,
    ).stdout.strip()
    (RESULTS_DIR / "booldog_commit.txt").write_text(booldog_commit + "\n")
    print(f"recorded biodivine commit ({commit}) and BoolDog commit "
          f"({booldog_commit}) in {RESULTS_DIR}")


def run_indegree_benchmark():
    items = [({"k": k, "in_degree": k + 2}, synthetic_indegree_bnet(k))
             for k in INDEGREE_VALUES]
    rows = run_sweep(items, stop_on_timeout=True)
    write_tsv(rows, ["k", "in_degree", *FIELDNAMES], RESULTS_DIR / "indegree.tsv")


def run_all():
    run_nodecount_benchmark()
    run_indegree_benchmark()


# --------------------------------------------------------------------
# plot
# --------------------------------------------------------------------

def load_tsv(name):
    with open(RESULTS_DIR / name) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return pd.DataFrame(rows)


def load_nodecount():
    df = load_tsv("nodecount.tsv")
    df = df[df["error"] == ""].copy()
    for col in ("n_nodes", "max_in_degree", *METRIC_LABELS):
        df[col] = pd.to_numeric(df[col])
    return df


def load_indegree():
    df = load_tsv("indegree.tsv")
    df["in_degree"] = pd.to_numeric(df["in_degree"])
    completed = df[df["error"] == ""].copy()
    for col in METRIC_LABELS:
        completed[col] = pd.to_numeric(completed[col])
    timed_out = df[df["error"] == "timeout"]
    first_timeout = timed_out.iloc[0]["in_degree"] if len(timed_out) else None
    return completed, first_timeout


def scatter_metrics(ax, x, color_col, df, cmap, norm, connect=False):
    df = df.sort_values(x)
    for col, label in METRIC_LABELS.items():
        if connect:
            ax.plot(df[x], df[col], color="0.6", linewidth=1, zorder=1)
        ax.scatter(
            df[x], df[col], c=df[color_col], cmap=cmap, norm=norm,
            marker=MARKERS_BY_METRIC[label], s=110, edgecolor="black", linewidth=0.5,
            zorder=2,
        )


def plot_combined():
    nodecount = load_nodecount()
    indegree, first_timeout = load_indegree()

    color_min = min(nodecount["max_in_degree"].min(), indegree["in_degree"].min())
    color_max = max(nodecount["max_in_degree"].max(), indegree["in_degree"].max())
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(color_min, color_max)

    y_min = min(nodecount[list(METRIC_LABELS)].min().min(),
                indegree[list(METRIC_LABELS)].min().min())
    y_max = max(nodecount[list(METRIC_LABELS)].max().max(),
                indegree[list(METRIC_LABELS)].max().max())

    plt.rcParams.update({"axes.grid": True, "grid.color": "0.85", "axes.facecolor": "white"})
    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    for ax in (ax_left, ax_right):
        ax.set_axisbelow(True)

    scatter_metrics(ax_left, "n_nodes", "max_in_degree", nodecount, cmap, norm)
    ax_left.set_xscale("log")
    ax_left.set_xlabel("Number of nodes")
    ax_left.set_title("Benchmark Dataset\n(biodivine-boolean-models)")

    scatter_metrics(ax_right, "in_degree", "in_degree", indegree, cmap, norm, connect=True)
    ax_right.set_xlabel("Max in-degree")
    ax_right.set_title("Controlled in-degree\n(synthetic models)")
    if first_timeout is not None:
        ax_right.axvline(first_timeout, color="firebrick", linestyle="--", alpha=0.7)
        ax_right.text(first_timeout, y_min, f"timeout at in-degree={first_timeout}  ",
                      color="firebrick", va="bottom", ha="right", fontsize=9, rotation=90)

    ax_left.set_yscale("log")
    ax_left.set_ylim(y_min / 2, y_max * 2)
    ax_left.set_ylabel("Runtime (s, log scale)")

    for _, label in METRIC_LABELS.items():
        ax_left.scatter([], [], facecolor="gray", edgecolor="black",
                        marker=MARKERS_BY_METRIC[label], s=110, label=label)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=ax_right, label="Max in-degree", pad=0.02)

    # Placed to the right of the colorbar
    metric_legend = ax_left.legend(
        title=None, bbox_to_anchor=(1.20, 0.5), bbox_transform=ax_right.transAxes,
        loc="center left",
    )

    # Commit-provenance box, placed directly below the marker legend
    biodivine_commit = (RESULTS_DIR / "biodivine_commit.txt").read_text().strip()
    booldog_commit = (RESULTS_DIR / "booldog_commit.txt").read_text().strip()
    commit_text = (
        f"biodivine-boolean-models @ {biodivine_commit[:10]}\n"
        f"BoolDog @ {booldog_commit[:10]}"
    )
    fig.canvas.draw()
    legend_bbox_fig = metric_legend.get_window_extent(fig.canvas.get_renderer()) \
        .transformed(fig.transFigure.inverted())
    legend_frame = metric_legend.get_frame()
    fontsize = metric_legend.get_texts()[0].get_fontsize()

    # Match the legend frame's own boxstyle exactly
    borderpad = plt.rcParams["legend.borderpad"]
    pad_fig_x = (borderpad * fontsize / 72.0) / fig.get_figwidth()
    fig.text(
        legend_bbox_fig.x0 + pad_fig_x, legend_bbox_fig.y0 - 0.03, commit_text,
        va="top", ha="left", fontsize=fontsize,
        family=metric_legend.get_texts()[0].get_fontfamily(),
        bbox=dict(boxstyle=f"round,pad={borderpad},rounding_size=0.2",
                  facecolor=legend_frame.get_facecolor(),
                  edgecolor=legend_frame.get_edgecolor(), linewidth=legend_frame.get_linewidth()),
    )

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES_DIR / "scaling.png", dpi=150, bbox_inches="tight",
                bbox_extra_artists=[metric_legend])
    print(f"wrote {FIGURES_DIR / 'scaling.png'}")


# --------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["run", "plot", "all"])
    args = parser.parse_args()

    if args.command in ("run", "all"):
        run_all()
    if args.command in ("plot", "all"):
        plot_combined()
