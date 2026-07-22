# Releasing BoolDog

Releases are made via the [`Release`](.github/workflows/release.yml)
GitHub Actions workflow, which chains into
[`python-publish.yml`](.github/workflows/python-publish.yml) to publish to
PyPI. Zenodo archiving (a distinct DOI per release) happens automatically
via Zenodo's own GitHub integration, triggered by the same GitHub Release.

## Steps

1. Go to the repo's **Actions** tab -> **Release** workflow -> **Run workflow**.
2. Choose the version bump (`patch`, `minor`, or `major`) and run it.
3. The workflow, in order:
   - runs the test suite (`test` job) -- the release is aborted if this fails;
   - bumps `version` in `pyproject.toml`, commits and pushes it to `main`;
   - creates and pushes a `X.Y.Z` tag;
   - creates a GitHub Release (name `X.Y.Z`, tag `X.Y.Z`, auto-generated
     release notes).
4. That Release being published triggers `python-publish.yml`, which
   builds the sdist/wheel and publishes them to PyPI via
   [trusted publishing](https://docs.pypi.org/trusted-publishers/)
   (no API token needed -- configured on PyPI's project settings page,
   under **Publishing**, matching this repo, `python-publish.yml`, and
   the `pypi` environment name).
5. Zenodo picks up the new GitHub Release independently and archives it
   under a new DOI (configured via Zenodo's GitHub integration, not
   anything in this repo).

Equivalent via the `gh` CLI instead of the Actions tab:

```bash
gh workflow run release.yml -f bump=patch
```

## Manual fallback

If the workflow can't be used for some reason, the equivalent manual
steps are:

```bash
# 1. bump `version` in pyproject.toml by hand, commit, push to main
# 2. tag and push
git tag X.Y.Z
git push origin X.Y.Z
# 3. create the GitHub release (triggers python-publish.yml)
gh release create X.Y.Z --title "X.Y.Z" --generate-notes
```

(Or, to publish to PyPI directly without going through GitHub Actions at
all: `python -m pip install --upgrade build twine`, `python -m build`,
`twine upload dist/*` -- requires a PyPI API token, since trusted
publishing only applies to the GitHub Actions workflow.)
