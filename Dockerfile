# BoolDog depends on pyboolnet, which bundles precompiled x86_64-only
# binaries (BNetToPrime, clasp, gringo, ...) with no ARM64 build. On an
# Apple Silicon host, pass `--platform linux/amd64` to `docker build` and
# `docker run` so Docker emulates x86_64 instead of picking a native arm64
# base image, for which pyboolnet has no matching binaries at all (see
# docs/source/installation.rst, "Known installation issues").
FROM python:3.12-slim

# cpp: eqntott_linux64 shells out to a hardcoded /lib/cpp; installing the
# cpp package creates that symlink automatically via update-alternatives.
# graphviz, imagemagick: used by BoolDog's STG/simulation-result plotting.
RUN apt-get update && apt-get install -y --no-install-recommends \
        cpp \
        graphviz \
        imagemagick \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .
RUN pip install --no-cache-dir .[all]

CMD ["python"]
