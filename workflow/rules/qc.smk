# ---------------------------------------------------------------------------
# qc.smk: input validation and quality control.
# ---------------------------------------------------------------------------


rule validate:
    """Check the inputs before anything expensive runs.

    Every check corresponds to a real silent failure from the analyses this
    method was extracted from. Fails the run rather than warning, because the
    failure mode these guard against is a confident wrong answer, not a crash.
    """
    input:
        samplesheet=config["input"]["samplesheet"],
        counts=config["input"]["counts"],
        **({"features": config["input"]["features"]} if REGION_MODE else {}),
    output:
        report=P("validation"),
        libsizes=P("libsizes"),
    log:
        f"{OUT}/logs/validate.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/01_validate.R"
