#!/usr/bin/env python3
"""Deep-merge YAML override files onto a base config and print the result.

    python tests/helpers/merge_config.py base.yaml override.yaml [...] > out.yaml

Concatenating YAML files is not a merge: a repeated top-level key replaces
the whole block, so an override of `input.features` silently discards
`input.counts` and `input.samplesheet`. That bug is why this exists.
"""
import sys
import yaml


def merge(a, b):
    for k, v in (b or {}).items():
        a[k] = merge(a.get(k) or {}, v) if isinstance(v, dict) and isinstance(a.get(k), dict) else v
    return a


cfg = yaml.safe_load(open(sys.argv[1])) or {}
for path in sys.argv[2:]:
    cfg = merge(cfg, yaml.safe_load(open(path)))
yaml.safe_dump(cfg, sys.stdout, sort_keys=False)
