'''Model definitions and fetching for the performance/scalability
benchmark: the curated list of real biological models (node-count axis,
downloaded on demand -- see fetch_all()), and a synthetic-model generator
(in-degree axis).
'''
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

BIODIVINE_REPO = "sybila/biodivine-boolean-models"
BIODIVINE_REF = "main"

DATA_DIR = Path(__file__).parent / "data"
BIODIVINE_COMMIT_FILE = DATA_DIR / "biodivine_commit.txt"


def max_in_degree(bnet_text):
    '''Maximum number of distinct regulators of any single rule in the
    given bnet text -- a cheap proxy computed from the raw text, so it's
    available even if prime-implicant conversion itself times out.
    '''
    lines = bnet_text.splitlines()
    if lines and lines[0].strip().lower().replace(" ", "") == "targets,factors":
        lines = lines[1:]
    max_arity = 0
    for line in lines:
        if "," not in line or line.strip().startswith("#"):
            continue
        _, rule = line.split(",", 1)
        literals = set(re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", rule))
        literals -= {"and", "or", "not"}
        max_arity = max(max_arity, len(literals))
    return max_arity


def synthetic_indegree_bnet(k):
    '''bnet text for a single synthetic node with in-degree k+2, reproducing
    the exact clause pattern from the known pyboolnet performance issue
    (hklarner/pyboolnet#122): an OR of (activator AND shared inhibitors)
    clauses, e.g. the PAX5 rule in that issue. Isolates the
    in-degree-vs-runtime relationship from network size.
    '''
    lines = [f"A{i}, A{i}" for i in range(k)]  # free/input nodes
    lines += ["M, M", "X, X"]
    clauses = " | ".join(f"(!M & A{i} & !X)" for i in range(k))
    lines.append(f"TARGET, {clauses}")
    return "\n".join(lines)


# label: short, unique identifier used in results/plots (biodivine's own
#   model id, stable and traceable back to its metadata.json/README).
# dirname: the model's directory name in the upstream repo's models/ folder.
#
# Stratified by node count (biodivine's own models/summary.csv, sampled
# roughly every 8th percentile) across the full size range available, PLUS
# a few additional models in the 90-320 node range deliberately chosen for
# LOW max in-degree (id-012, id-016, id-221) alongside the size-stratified
# picks that happen to have high in-degree in the same range (id-087,
# id-151, id-239) -- same node-count neighbourhood, very different in-degree
# and runtime, which is the whole point of the hue-coloured scaling plot.
MODELS = [
    {"label": "id-165", "dirname": "[id-165]__[var-4]__[in-4]__[EGGSHELL-PATTERNING-PHENOMOENOLOGICAL]"},
    {"label": "id-095", "dirname": "[id-095]__[var-9]__[in-1]__[FISSION-YEAST-2008]"},
    {"label": "id-281", "dirname": "[id-281]__[var-12]__[in-0]__[EMT-SWITCH]"},
    {"label": "id-164", "dirname": "[id-164]__[var-17]__[in-7]__[EGGSHELL-PATTERNING-MECHANISTIC]"},
    {"label": "id-069", "dirname": "[id-069]__[var-20]__[in-2]__[IRON-ACQUISITION-AND-STRESS-RESPONSE]"},
    {"label": "id-005", "dirname": "[id-005]__[var-28]__[in-0]__[FA-BRCA-PATHWAY]"},
    {"label": "id-129", "dirname": "[id-129]__[var-33]__[in-1]__[RTC-AND-TRANSCRIPTION]"},
    {"label": "id-060", "dirname": "[id-060]__[var-44]__[in-5]__[STOMATAL-OPENING]"},
    {"label": "id-259", "dirname": "[id-259]__[var-53]__[in-23]__[IMMUNE-SURVEILLANCE-OF-PRIMARY-MELANOMA-LC]"},
    {"label": "id-226", "dirname": "[id-226]__[var-68]__[in-7]__[B-CELL-APOPTOSIS]"},
    {"label": "id-012", "dirname": "[id-012]__[var-94]__[in-7]__[T-CELL-RECEPTOR-SIGNALING]"},
    {"label": "id-087", "dirname": "[id-087]__[var-99]__[in-34]__[INFLAMMATORY-GENE-EXPRESSION-IN-MACROPHAGES]"},
    {"label": "id-016", "dirname": "[id-016]__[var-104]__[in-14]__[IL-1-SIGNALING]"},
    {"label": "id-151", "dirname": "[id-151]__[var-130]__[in-3]__[TCR-REDOX-METABOLISM]"},
    {"label": "id-239", "dirname": "[id-239]__[var-234]__[in-69]__[M1-SYNOVIAL-MACROPHAGE]"},
    {"label": "id-221", "dirname": "[id-221]__[var-280]__[in-37]__[MYCOBACTERIAL-LATENCY]"},
    {"label": "id-001", "dirname": "[id-001]__[var-302]__[in-19]__[SIGNALING-IN-MACROPHAGE-ACTIVATION]"},
    {"label": "id-252", "dirname": "[id-252]__[var-695]__[in-65]__[MAMMALIAN-EPIDERMIS-2D]"},
    {"label": "id-243", "dirname": "[id-243]__[var-853]__[in-223]__[RHEUMATOID-ARTHRITIS-MULTI-CELLULAR]"},
]


def resolve_commit():
    '''Resolve BIODIVINE_REF to a concrete commit SHA via the GitHub API,
    so every model file in a given fetch_all() call comes from the same
    snapshot even if upstream changes mid-fetch.
    '''
    url = f"https://api.github.com/repos/{BIODIVINE_REPO}/commits/{BIODIVINE_REF}"
    with urllib.request.urlopen(url) as response:
        return json.load(response)["sha"]


def fetch_all():
    '''Download the bnet files listed in MODELS into data/ (not committed
    to the BoolDog repository, see .gitignore).

    Fetches from the current BIODIVINE_REF branch rather than a pinned
    commit, but resolves and records the exact commit SHA actually used
    (data/biodivine_commit.txt), so a benchmark run's report can cite
    precisely what it fetched even though the fetch itself always tracks
    upstream.
    '''
    DATA_DIR.mkdir(exist_ok=True)

    if BIODIVINE_COMMIT_FILE.exists():
        commit = BIODIVINE_COMMIT_FILE.read_text().strip()
        print(f"using previously recorded commit: {commit}")
    else:
        commit = resolve_commit()
        BIODIVINE_COMMIT_FILE.write_text(commit + "\n")
        print(f"resolved {BIODIVINE_REPO}@{BIODIVINE_REF} -> {commit}")

    for model in MODELS:
        dest = DATA_DIR / f"{model['label']}.bnet"
        if dest.exists():
            print(f"skip (cached): {dest.name}")
            continue
        url = (
            f"https://raw.githubusercontent.com/{BIODIVINE_REPO}/"
            f"{commit}/models/{urllib.parse.quote(model['dirname'])}/model.bnet"
        )
        print(f"fetching {model['label']} <- {url}")
        with urllib.request.urlopen(url) as response:
            dest.write_bytes(response.read())

    return commit


if __name__ == "__main__":
    fetch_all()
