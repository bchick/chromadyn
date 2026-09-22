#!/usr/bin/env python3
"""Every cfg("a.b.c") key used by a script must exist in config.yaml.

A key that does not exist silently falls back to the default baked into the
call, which is a threshold living in a script: exactly what this project
forbids.
"""
import pathlib
import re
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


known = keys(yaml.safe_load(open("config/config.yaml")))
pat = re.compile(r'\bcfg(?:_req)?\(\s*"([a-z0-9_.]+)"')
ok = True
for f in sorted(pathlib.Path("workflow/scripts").rglob("*.R")):
    for m in pat.finditer(f.read_text()):
        if m.group(1) not in known:
            print(f"{f}: cfg(\"{m.group(1)}\") is not a key in config/config.yaml", file=sys.stderr)
            ok = False
sys.exit(0 if ok else 1)
