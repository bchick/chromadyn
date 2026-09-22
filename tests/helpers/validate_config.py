#!/usr/bin/env python3
"""Validate one config file against config/schema.yaml, merged over the default.

A preset or demo config is a partial override, so it is validated as the
merged result rather than on its own.
"""
import sys
import jsonschema
import yaml


def merge(a, b):
    for k, v in b.items():
        a[k] = merge(a.get(k, {}), v) if isinstance(v, dict) and isinstance(a.get(k), dict) else v
    return a


def main(path):
    base = yaml.safe_load(open("config/config.yaml"))
    cfg = base if path == "config/config.yaml" else merge(base, yaml.safe_load(open(path)) or {})
    jsonschema.validate(cfg, yaml.safe_load(open("config/schema.yaml")))


if __name__ == "__main__":
    main(sys.argv[1])
