# CoLoMoTo Docker integration (exploratory)

Stage 0 of the [CoLoMoTo integration
process](https://github.com/colomoto/colomoto-docker/blob/for-next/CONTRIBUTING.md):
an extension of the existing
[colomoto-docker](https://github.com/colomoto/colomoto-docker) image.

## Usage

Build the image from the repo root:

```bash
docker build -t booldog-colomoto colomoto/
```

Then start a notebook server with the repo bound in, e.g. via the
[`colomoto-docker` launcher](https://pypi.org/project/colomoto-docker/):

```bash
colomoto-docker --image booldog-colomoto --no-update --bind .
```

or with plain `docker run` (passing your host UID/GID so bind-mounted files
keep your ownership, not root's):

```bash
docker run --rm -p 8888:8888 \
  --volume "$(pwd)":/notebook -w /notebook \
  -e UID=$(id -u) -e GID=$(id -g) \
  booldog-colomoto \
  jupyter-notebook --no-browser --port 8888 --ip 0.0.0.0 --NotebookApp.token=
```

Open `http://127.0.0.1:8888/tree` and navigate to `tutorials/` for BoolDoG's
tutorial notebooks.

## Notes

⚠️ Unlike CONTRIBUTING.md's own Stage 0 template, this Dockerfile does *not*
end with `USER user`. The base image's entrypoint needs to start as root so
it can `gosu user ...` to drop to the UID/GID passed in above; a trailing
`USER user` in the Dockerfile pre-empts that and breaks it with
`error: failed switching to 'user': operation not permitted`.

BoolDog works correctly in this image, using the same `pyboolnet` version
already bundled there.

⚠️ Installing the `sbml` extra (`booldog[all]` or `booldog[sbml]`) breaks
**NORDic** and **AstroLogics**, two tools already in this image. If
you need those tools alongside BoolDog, install without `sbml`:

```bash
pip install booldog[networks,graphviz,biomodels]
```

The cause is upstream SBML specific BoolDog dependencies (via `tabularqual`/`sbmlutils`/`py4cytoscape`) forcing a `numpy` upgrade (`2.2.0` → `2.5.1`). Solution proposed at [matthiaskoenig/sbmlutils#461](https://github.com/matthiaskoenig/sbmlutils/issues/461).
