# ---------------------------------------------------------------------------
# differential.smk: per-arm modelling and dynamic-feature selection.
#
# Both rules fan out over the `arm` wildcard. Each arm is fitted separately
# against the baseline libraries it shares with the others, which is why
# arm_samples() rather than a group filter decides the columns.
# ---------------------------------------------------------------------------


rule differential:
    """Fit one model per arm: LRT over time, per-timepoint Wald contrasts, or none."""
    input:
        dds=P("obj_dds"),
        transformed=P("obj_transformed"),
    output:
        obj=P("obj_differential"),
    params:
        samples=lambda w: arm_samples(w.arm),
    log:
        f"{OUT}/logs/differential_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/05_differential.R"


rule select_dynamic:
    """Apply both gates: significance, and a floor on the range of timepoint means."""
    input:
        obj=P("obj_differential"),
    output:
        results=P("differential"),
    log:
        f"{OUT}/logs/select_dynamic_{{arm}}.log",
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/06_select_dynamic.R"
