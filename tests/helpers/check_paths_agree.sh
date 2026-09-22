#!/usr/bin/env bash
# Resolve every workflow/paths.yaml template on the Python side and the R side
# and assert the two agree. Both read the same file, but through different
# formatters, and a divergence would not surface until a rule wrote its output
# somewhere the next rule does not look.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1

pyout=$(python - <<'PY'
import yaml
p = yaml.safe_load(open("workflow/paths.yaml"))
for k in sorted(p):
    print(k, p[k].format(out="results", arm="ARM", fmt="png", sc="SC",
                         rule="RULE", arm_suffix="_ARM"))
PY
) || exit 1

rout=$(Rscript - <<'RS'
p <- yaml::read_yaml("workflow/paths.yaml")
subs <- c(out="results", arm="ARM", fmt="png", sc="SC", rule="RULE", arm_suffix="_ARM")
for (k in sort(names(p))) {
  t <- p[[k]]
  for (n in names(subs)) t <- gsub(paste0("{", n, "}"), subs[[n]], t, fixed = TRUE)
  cat(k, t, "\n")
}
RS
) || exit 1

diff <(echo "$pyout") <(echo "$rout" | sed 's/ $//') >&2
