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


rule prefilter:
    """Build the DESeqDataSet and drop features below filter.min_mean_counts."""
    input:
        samplesheet=config["input"]["samplesheet"],
        counts=config["input"]["counts"],
        validation=P("validation"),
    output:
        dds=P("obj_dds"),
    log:
        f"{OUT}/logs/prefilter.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/02_prefilter.R"


rule transform:
    """Variance-stabilize, and optionally remove a batch effect.

    The output of this rule is the common scale everything downstream works
    on, which is why differential.min_range is expressed against it.
    """
    input:
        dds=P("obj_dds"),
    output:
        transformed=P("obj_transformed"),
    log:
        f"{OUT}/logs/transform.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/03_transform.R"


rule qc:
    """PCA, sample correlation and replicate correlation.

    Advisory: nothing here gates the run. It exists so a design problem is
    visible before the clusters built on it get interpreted.
    """
    input:
        dds=P("obj_dds"),
        transformed=P("obj_transformed"),
    output:
        rep_cor_table=P("replicate_cor_tsv"),
        figures=expand(P("pca"), fmt=config["figures"]["formats"])
               + expand(P("sample_cor"), fmt=config["figures"]["formats"])
               + expand(P("replicate_cor"), fmt=config["figures"]["formats"]),
    log:
        f"{OUT}/logs/qc.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/04_qc.R"
