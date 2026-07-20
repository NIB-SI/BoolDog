# Stage 0 of the CoLoMoTo integration process (see CONTRIBUTING.md at
# https://github.com/colomoto/colomoto-docker/blob/for-next/CONTRIBUTING.md):
# a throwaway extension of the existing colomoto-docker image, used to test
# that BoolDog installs and runs correctly alongside the tools already
# bundled there (notably pyboolnet==3.0.16, which BoolDog also depends on),
# before attempting the real integration (a colomoto-conda recipe, a demo
# notebook, and a PR to colomoto-docker itself).
FROM colomoto/colomoto-docker:next

# Note: no trailing `USER user` here (unlike CONTRIBUTING.md's own template).
# The base image's entrypoint needs to start as root so it can `gosu user ...`
# to drop privileges matching the host UID/GID passed in by the colomoto-docker
# launcher's --bind; pre-emptively switching USER breaks that with
# "error: failed switching to 'user': operation not permitted".
USER root
RUN pip install --no-cache-dir booldog[all]
