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


rule superclusters:
    """Collapse clusters into named trajectory classes, and split one of them.

    Emits the dendrogram fit and a k-selection diagnostic whatever the
    assignment method is, so the cut that was taken can be checked against
    the cut that was not.
    """
    input:
        obj=P("obj_degpatterns"),
        **({"features": config["input"]["features"]} if REGION_MODE else {}),
    output:
        clusters=P("clusters"),
        profiles=P("cluster_profiles"),
        sizes=P("supercluster_sizes"),
        kdiag=P("kdiag"),
        fit=P("obj_supercluster"),
    log:
        f"{OUT}/logs/superclusters_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/08_superclusters.R"


rule figures:
    """Every per-arm figure.

    Declares a single flag file rather than every format of every panel, so
    that changing figures.formats does not multiply the DAG. The flag lists
    what was written.
    """
    input:
        cluster_obj=P("obj_degpatterns"),
        fit=P("obj_supercluster"),
        results=P("differential"),
        transformed=P("obj_transformed"),
        dds=P("obj_dds"),
    output:
        flag=P("fig_flag"),
    log:
        f"{OUT}/logs/figures_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/10_figures.R"
