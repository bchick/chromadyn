# ---------------------------------------------------------------------------
# chromadyn container.
#
#   docker build -t chromadyn .
#   docker run --rm -v "$PWD:/work" -w /work chromadyn \
#       snakemake --configfile config/demo.yaml -j 4
#
# Built on the pixi image so the environment is exactly what pixi.lock pins:
# the container, `pixi install` and `snakemake --use-conda` are three routes
# to one stack rather than three stacks that drift.
# ---------------------------------------------------------------------------
FROM ghcr.io/prefix-dev/pixi:0.64.0 AS build

WORKDIR /opt/chromadyn

# Resolve the environment first, so a change to the workflow does not
# invalidate the dependency layer.
COPY pixi.toml pixi.lock ./
RUN pixi install --locked

COPY workflow/ ./workflow/
COPY config/ ./config/
COPY demo/ ./demo/
COPY tests/ ./tests/
COPY LICENSE README.md CITATION.cff ./

# Complete the Bioconductor install that pixi leaves undone; see
# workflow/envs/postinstall.sh for why this is a separate step.
RUN pixi run postinstall

# `pixi shell-hook` bakes the activation into a plain entrypoint, so the
# container does not need pixi at run time.
RUN pixi shell-hook -e default > /shell-hook.sh && \
    echo 'exec "$@"' >> /shell-hook.sh

FROM ghcr.io/prefix-dev/pixi:0.64.0 AS runtime
WORKDIR /opt/chromadyn
COPY --from=build /opt/chromadyn /opt/chromadyn
COPY --from=build /shell-hook.sh /shell-hook.sh

LABEL org.opencontainers.image.title="chromadyn" \
      org.opencontainers.image.description="Temporal clustering of omics timecourses" \
      org.opencontainers.image.source="https://github.com/bchick/chromadyn" \
      org.opencontainers.image.licenses="MIT"

ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]
CMD ["snakemake", "--configfile", "config/demo.yaml", "-j", "4"]
