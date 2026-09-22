# Contributing

## The rules that matter

**Every threshold lives in the config.** No numeric threshold may appear in a
script. This is the most-violated rule in the project this method came from
and the main reason it was not reusable. `tests/test_static.sh` greps for
violations and fails the build. If you need a new threshold, add it to
`config/config.yaml` with a comment explaining *why the key exists*, add it to
`config/schema.yaml`, and read it with `cfg()`.

**Every output path lives in `workflow/paths.yaml`.** Both the Snakefile and
the R library read it, and tier 1 diffs the two resolutions. Never build a
path by pasting strings.

**No em dashes** anywhere in prose, docs, comments or commit messages. Commas,
parentheses, semicolons, colons, or split the sentence. En dashes in numeric
ranges are fine. Tier 1 enforces this.

**No AI attribution** in commits, pull requests, README bylines, `CITATION.cff`
authors or code headers. This work is published alongside academic
manuscripts. Tier 1 enforces this too.

## R style

Follow the surrounding code. `.lintr` holds the enforced subset, which is the
correctness linters rather than the cosmetic ones; run `pixi run --frozen
Rscript tests/helpers/run_lintr.R`.

Specifics worth stating:

- Every script opens with the same header: the `if (!exists("snakemake"))`
  standalone block, then `lib/common.R`, then `lib/namespace_fix.R`, then
  `cd_init(snakemake)`. Keep them identical, so a reader who has read one has
  read all of them.
- The standalone block names `workflow/paths.yaml` keys, not literal paths, so
  a debug entry point cannot drift from the rule it mirrors.
- `lib/namespace_fix.R` must be sourced non-locally and must assign into
  `.GlobalEnv`. Sourcing it inside a function does nothing.
- Palettes live in `lib/theme.R` and nowhere else. The source project
  redefined one in a notebook after sourcing the shared file, with different
  hex values, and two panels of the same figure disagreed.
- Threshold and reference lines use `col_accent` (`#e6ab02`), never red.
- `save_figure()` takes **pixels at 72 per inch**, not inches. Use the
  `FIG_W_*` constants.

## Shell style

`set -uo pipefail`, deliberately without `-e`, so a failing assertion records
and the suite continues. Source `tests/lib/assert.sh` for the
`pass/fail/warn/skip/group` helpers. `shellcheck -S warning` must be clean.

## Before opening a pull request

```bash
bash tests/test_static.sh          # seconds
bash tests/test_endtoend.sh        # minutes
Rscript tests/test_correctness.R   # minutes
```

If you changed the method deliberately, the golden tables will differ.
Regenerate them with `UPDATE_GOLDEN=1` and **say what changed and why in the
pull request**. A golden diff that is not explained is a regression.

## Adding a rule

1. Add its output templates to `workflow/paths.yaml`.
2. Add the rule to the right `workflow/rules/*.smk`, with `conda:` and `log:`.
3. Gate optional rules in `final_targets()`, never by wrapping the `rule:`
   block in an `if`. The linter and tier 1 must be able to see every rule.
4. Write the script with the standard header.
5. Add its outputs to tier 2's existence checks.

## Commits and releases

Conventional prefixes: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`,
`chore:`. Write commit messages that say why, not what; the diff already says
what.

Releases are semver tags. Tagging `v*` builds and pushes the container to
GHCR. For a citable DOI, enable the Zenodo GitHub integration for this
repository before tagging: Zenodo archives each release and mints a DOI, and
the concept DOI should then be added to `CITATION.cff`.
