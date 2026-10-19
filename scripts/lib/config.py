#!/usr/bin/env python3
"""Read scalar values from the repository YAML configuration."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
from typing import Any


def nested_get(value: Any, dotted_key: str) -> Any:
    for part in dotted_key.split("."):
        if not isinstance(value, dict) or part not in value:
            raise KeyError(dotted_key)
        value = value[part]
    return value


def parse_scalar(raw: str) -> Any:
    value = raw.strip()
    if not value:
        return {}
    if value in {"null", "Null", "NULL", "~"}:
        return None
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    try:
        return ast.literal_eval(value)
    except (ValueError, SyntaxError):
        return value


def load_nested_scalars(path: Path) -> dict[str, Any]:
    """Parse the indented mapping subset used by this repository.

    Lists are intentionally not supported because shell wrappers query only scalar
    configuration values. Full YAML parsing remains available to the R modules.
    """
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.lstrip(" ")
        if not stripped or stripped.startswith("#") or stripped.startswith("-"):
            continue
        indent = len(line) - len(stripped)
        if ":" not in stripped:
            raise ValueError(f"Unsupported YAML at line {number}")
        key, raw = stripped.split(":", 1)
        while stack[-1][0] >= indent:
            stack.pop()
        parent = stack[-1][1]
        parsed = parse_scalar(raw)
        parent[key.strip()] = parsed
        if isinstance(parsed, dict):
            stack.append((indent, parsed))
    return root


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("key")
    parser.add_argument("--path", action="store_true", help="Resolve a repository-relative path")
    args = parser.parse_args()

    config_path = args.config.resolve()
    config = load_nested_scalars(config_path)
    value = nested_get(config, args.key)
    if isinstance(value, (dict, list)):
        raise SystemExit(f"Configuration value is not scalar: {args.key}")
    if value is None:
        raise SystemExit(f"Configuration value is null: {args.key}")
    if args.path:
        path = Path(str(value))
        if not path.is_absolute():
            path = config_path.parent.parent / path
        value = path.resolve()
    if isinstance(value, bool):
        print(str(value).lower())
    else:
        print(value)


if __name__ == "__main__":
    main()
