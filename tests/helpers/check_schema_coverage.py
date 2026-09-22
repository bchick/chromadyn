#!/usr/bin/env python3
"""Assert config.yaml and schema.yaml describe the same set of keys.

Catches documentation drift in both directions: a key the schema requires but
the shipped config never shows (so nobody knows it exists), and a key the
config ships that the schema does not constrain (so a typo in it is silent).
"""
import sys
import yaml


def keys(node, prefix=""):
    out = set()
    if isinstance(node, dict):
        for k, v in node.items():
            p = f"{prefix}.{k}" if prefix else k
            out.add(p)
            out |= keys(v, p)
    return out


def schema_keys(node, prefix=""):
    out = set()
    props = node.get("properties", {})
    for k, v in props.items():
        p = f"{prefix}.{k}" if prefix else k
        out.add(p)
        out |= schema_keys(v, p)
        for branch in v.get("anyOf", []):
            out |= schema_keys(branch, p)
    return out


cfg = yaml.safe_load(open("config/config.yaml"))
sch = yaml.safe_load(open("config/schema.yaml"))

# Only compare the structural keys, not the contents of free-form maps such as
# superclusters.overrides, whose keys are user data rather than schema.
FREEFORM = {"superclusters.overrides"}
ck = {k for k in keys(cfg) if not any(k.startswith(f + ".") for f in FREEFORM)}
sk = schema_keys(sch)

undocumented = sk - ck          # schema knows it, config never shows it
unvalidated = ck - sk           # config ships it, schema ignores it

ok = True
for k in sorted(undocumented):
    print(f"schema declares '{k}' but config/config.yaml does not ship it", file=sys.stderr)
    ok = False
for k in sorted(unvalidated):
    print(f"config/config.yaml ships '{k}' but the schema does not constrain it", file=sys.stderr)
    ok = False
sys.exit(0 if ok else 1)
