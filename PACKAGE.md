# Releasing BoolDog

Releases are made via the single [`Release`](.github/workflows/release-and-publish.yml)
GitHub Actions workflow -- test, version bump, tag, GitHub Release, build,
and PyPI publish are all jobs in one `workflow_dispatch` run, chained via
`needs:`.

Zenodo archiving (a distinct DOI per release) happens automatically via
Zenodo's own GitHub integration, triggered by the same GitHub Release.

## Steps

1. Go to the repo's **Actions** tab -> **Release** workflow -> **Run workflow**.
2. Choose the version bump (`patch`, `minor`, or `major`) and run it.
3. In order: `test` (aborts the release if it fails) -> `release` (bumps
   `version` in `pyproject.toml`, commits and pushes to `main`, tags
   `X.Y.Z`, creates the GitHub Release) -> `build` (sdist + wheel) ->
   `pypi-publish` (publishes via
   [trusted publishing](https://docs.pypi.org/trusted-publishers/), no
   API token needed -- configured on PyPI's project settings page, under
   **Publishing**, matching this repo, the workflow filename, and the
   `pypi` environment name).

Equivalent via the `gh` CLI instead of the Actions tab:

```bash
gh workflow run release-and-publish.yml -f bump=patch
```

## Manual fallback

If the workflow can't be used for some reason, the equivalent manual
steps are:

```bash
# 1. bump `version` in pyproject.toml by hand, commit, push to main
# 2. tag and push
git tag X.Y.Z
git push origin X.Y.Z
# 3. create the GitHub release
gh release create X.Y.Z --title "X.Y.Z" --generate-notes
# 4. build and publish to PyPI (trusted publishing only applies inside
#    the GitHub Actions workflow, so this needs a PyPI API token instead)
python -m pip install --upgrade build twine
python -m build
twine upload dist/*
```
