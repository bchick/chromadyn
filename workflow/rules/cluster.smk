# ---------------------------------------------------------------------------
# cluster.smk: trajectory clustering and trajectory-class assignment.
# ---------------------------------------------------------------------------


rule cluster:
    """Group dynamic features by trajectory shape, one run per arm.

    Caches the full DEGreport object so that re-running superclusters,
    annotation or figures never re-enters degPatterns, which is the expensive
    step by a wide margin.
    """
    input:
        dds=P("obj_dds"),
        transformed=P("obj_transformed"),
        diff_obj=P("obj_differential"),
        results=P("differential"),
    output:
        obj=P("obj_degpatterns"),
        qc=P("assignment_qc"),
    log:
        f"{OUT}/logs/cluster_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/07_cluster.R"
