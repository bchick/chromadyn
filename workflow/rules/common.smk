# ---------------------------------------------------------------------------
# common.smk: helpers shared by every rule file.
#
# Deliberately contains no rules. Everything here is pure: it reads the
# samplesheet and the path table and derives the values the DAG is built from,
# so that `snakemake -n` can show the whole graph without running anything.
# ---------------------------------------------------------------------------


def load_paths():
    """Output path templates, the same table workflow/scripts/lib/common.R reads."""
    with open(os.path.join(workflow.basedir, "paths.yaml")) as fh:
        return yaml.safe_load(fh)


def P(key, **kw):
    """Resolve a path template from workflow/paths.yaml.

    Leaves any placeholder not supplied here in place, so that Snakemake can
    treat it as a wildcard: P("clusters") -> "results/clusters/{arm}_clusters.tsv".
    """
    tmpl = PATHS.get(key)
    if tmpl is None:
        raise WorkflowError(
            f"'{key}' is not in workflow/paths.yaml. "
            f"Known keys: {', '.join(sorted(PATHS))}"
        )
    passthrough = {k: f"{{{k}}}" for k in _placeholders(tmpl)}
    return tmpl.format(out=OUT, **(passthrough | kw))


def _placeholders(tmpl):
    return {m for m in re.findall(r"\{(\w+)\}", tmpl) if m != "out"}


def load_samplesheet(path):
    """Read the samplesheet as a list of dicts, failing loudly on the basics.

    Full validation is 01_validate.R's job. This only reads enough to build the
    DAG, and only complains about things that would make the DAG itself wrong.
    """
    if not os.path.exists(path):
        raise WorkflowError(
            f"samplesheet not found: {path}\n"
            f"Set input.samplesheet in your config, or run from the repository root."
        )
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise WorkflowError(f"samplesheet is empty: {path}")

    required = ["sample", "group", "time", "replicate"]
    missing = [c for c in required if c not in rows[0]]
    if missing:
        raise WorkflowError(
            f"samplesheet {path} is missing column(s): {', '.join(missing)}. "
            f"Found: {', '.join(rows[0])}"
        )
    for i, r in enumerate(rows, start=2):
        try:
            float(r["time"])
        except (TypeError, ValueError):
            raise WorkflowError(
                f"samplesheet {path} line {i}: time '{r['time']}' is not numeric. "
                f"Units belong in config input.time_unit, never in the value."
            )
    return rows


def derive_arms(samples, baseline):
    """Treatment arms, in first-seen order. The baseline sentinel is not an arm.

    Rows whose group equals the sentinel join every arm rather than forming one
    of their own: that is how a multi-arm timecourse sharing one unstimulated
    baseline is expressed.
    """
    arms = list(dict.fromkeys(r["group"] for r in samples if r["group"] != baseline))
    if not arms:
        raise WorkflowError(
            f"no treatment arms: every samplesheet row has group == "
            f"input.baseline_group ('{baseline}'). Set baseline_group to a value "
            f"that only your baseline libraries use, or to something unused if "
            f"every arm carries its own baseline."
        )
    bad = [a for a in arms if not re.fullmatch(r"[A-Za-z0-9._+-]+", a)]
    if bad:
        raise WorkflowError(
            f"arm name(s) cannot be used as filename components: {', '.join(bad)}. "
            f"Use letters, digits, dot, underscore, plus or hyphen."
        )
    return arms


def arm_samples(arm):
    """Sample IDs for one arm, baseline included, ordered by time then replicate.

    Column order for every per-arm matrix is defined here and nowhere else, so
    that metadata built against it cannot drift out of alignment.
    """
    rows = [r for r in SAMPLES if r["group"] in (arm, BASELINE)]
    rows.sort(key=lambda r: (float(r["time"]), r["replicate"]))
    return [r["sample"] for r in rows]


def check_supercluster_cut(cfg):
    """k and h are alternative ways to cut one tree; both set is a contradiction.

    The schema rejects this too, but jsonschema's message for a failed `not` is
    close to unreadable, so say it properly here.
    """
    sc = cfg["superclusters"]
    if sc.get("k") is not None and sc.get("h") is not None:
        raise WorkflowError(
            f"superclusters.k ({sc['k']}) and superclusters.h ({sc['h']}) are both "
            f"set. They are alternative ways to cut the same tree: k gives a number "
            f"of classes, h a dendrogram height. Set exactly one and leave the other null."
        )
    if sc["method"] == "hclust" and sc.get("k") is None and sc.get("h") is None:
        raise WorkflowError(
            "superclusters.method is 'hclust' but neither k nor h is set, so there is "
            "no cut to make. Set one of them, or use method: shape."
        )
    if sc["method"] == "overrides" and not sc.get("overrides"):
        raise WorkflowError(
            "superclusters.method is 'overrides' but superclusters.overrides is empty. "
            "Supply a {cluster_id: name} map."
        )


def check_overrides(cfg):
    """Normalize the override map's keys to strings and check its values.

    YAML parses `{14: Decreasing}` with integer keys, which JSON Schema cannot
    express as property names, so the real check lives here. Whether the map
    covers every observed cluster is checked in R, where the clusters exist.
    """
    ov = cfg["superclusters"].get("overrides")
    if not ov:
        return
    bad = [k for k, v in ov.items() if not isinstance(v, str) or not v.strip()]
    if bad:
        raise WorkflowError(
            f"superclusters.overrides values must be non-empty class names; "
            f"offending key(s): {', '.join(map(str, bad))}"
        )
    cfg["superclusters"]["overrides"] = {str(k): v for k, v in ov.items()}


def final_targets(_):
    """Everything a complete run produces, given the configured mode.

    Optional rules are gated here, by whether their outputs are requested,
    rather than by wrapping the rule definitions in `if` blocks. That keeps
    every rule visible to the linter and to the tier 1 static test.
    """
    targets = [P("validation"), P("libsizes"), P("replicate_cor_tsv")]
    for key in ("pca", "sample_cor", "replicate_cor"):
        targets += expand(P(key), fmt=config["figures"]["formats"])
    targets += expand(P("differential"), arm=ARMS)
    targets += expand(P("assignment_qc"), arm=ARMS)
    for key in ("clusters", "cluster_profiles", "supercluster_sizes", "kdiag"):
        targets += expand(P(key), arm=ARMS)
    targets += expand(P("fig_flag"), arm=ARMS)
    return targets
