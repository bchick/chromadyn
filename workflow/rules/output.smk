# ---------------------------------------------------------------------------
# output.smk: exports, cross-arm comparison, report and manifest.
#
# export_bed and annotate are always DEFINED so that the linter and the tier 1
# static test can see them; whether they run is decided by final_targets(),
# not by wrapping these blocks in an `if`.
# ---------------------------------------------------------------------------


rule export_bed:
    """One BED6 per trajectory class. Region mode only.

    Uses directory() plus a manifest rather than a checkpoint: class names are
    data-dependent and contain spaces, so their filenames cannot be known at
    parse time, and nothing downstream consumes an individual BED.
    """
    input:
        clusters=P("clusters"),
    output:
        dir=directory(P("bed_dir")),
        manifest=P("bed_manifest"),
    log:
        f"{OUT}/logs/export_bed_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/09_export_bed.R"


rule annotate:
    """Nearest gene and distance to TSS. Needs annotate.txdb to be set."""
    input:
        clusters=P("clusters"),
    output:
        annot=P("annot"),
    log:
        f"{OUT}/logs/annotate_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/09b_annotate.R"


rule compare_arms:
    """Cross-tabulate class membership between arms.

    Depends on the cluster tables only. It deliberately does NOT depend on
    annotate: in gene mode that rule never runs, and making the comparison
    wait on it would serialize the whole tail of the workflow behind a rule
    that is not there.
    """
    input:
        clusters=expand(P("clusters"), arm=ARMS),
    output:
        crosstab=P("crosstab"),
    params:
        arms=ARMS,
    log:
        f"{OUT}/logs/compare_arms.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/11_compare_arms.R"


rule report:
    """Render the HTML report."""
    input:
        validation=P("validation"),
        libsizes=P("libsizes"),
        crosstab=P("crosstab"),
        clusters=expand(P("clusters"), arm=ARMS),
        sizes=expand(P("supercluster_sizes"), arm=ARMS),
        profiles=expand(P("cluster_profiles"), arm=ARMS),
        kdiag=expand(P("kdiag"), arm=ARMS),
        qc=expand(P("assignment_qc"), arm=ARMS),
        selection=expand(P("cluster_selection"), arm=ARMS),
        figures=expand(P("fig_flag"), arm=ARMS),
    output:
        html=P("report"),
    params:
        arms=ARMS,
    log:
        f"{OUT}/logs/report.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/13_report.R"


rule manifest:
    """Resolved config, versions, seeds and input checksums.

    Runs last: it merges the provenance fragment every other rule wrote as it
    finished, so it must see all of them.
    """
    input:
        report=P("report"),
    output:
        manifest=P("manifest"),
    log:
        f"{OUT}/logs/manifest.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/12_manifest.R"
