# coding: utf-8
"""定点修改 CNVkit params.py 中的六个硬编码常量。

不重写整个文件，也不修改 reference.py。
CNVkit 0.9.14 在 _reconcile_sex_guesses() 中比较 target 与 antitarget 的
chrX 证据。旧的 reference.py 字符串替换不再适用。
"""

from __future__ import annotations

import argparse
import ast
import difflib
import re
import shutil
import sys
from pathlib import Path

PARAM_TYPES = {
    "MIN_REF_COVERAGE": float,
    "MAX_REF_SPREAD": float,
    "NULL_LOG2_COVERAGE": float,
    "GC_MIN_FRACTION": float,
    "GC_MAX_FRACTION": float,
    "INSERT_SIZE": int,
}

ASSIGN_RE = re.compile(
    r"^(?P<name>" + "|".join(PARAM_TYPES) + r")\s*=\s*(?P<value>[^#\n]+)",
    re.MULTILINE,
)


def str2bool(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected a boolean, got %r" % (value,))


def default_params_path():
    import cnvlib

    return Path(cnvlib.__file__).resolve().parent / "params.py"


def read_assignments(text):
    found = {}
    for match in ASSIGN_RE.finditer(text):
        name = match.group("name")
        if name in found:
            raise SystemExit("params.py has repeated assignments: %s" % name)
        found[name] = match.group("value").strip()
    missing = [name for name in PARAM_TYPES if name not in found]
    if missing:
        raise SystemExit("params.py is missing assignments: %s" % ", ".join(missing))
    return found


def parse_literal(name, raw):
    try:
        value = ast.literal_eval(raw)
    except (SyntaxError, ValueError) as exc:
        raise SystemExit("cannot parse %s = %r: %s" % (name, raw, exc))
    expected = PARAM_TYPES[name]
    if isinstance(value, bool) or not isinstance(value, expected):
        raise SystemExit(
            "%s must be %s, got %s" % (name, expected.__name__, type(value).__name__)
        )
    return value


def validate(values):
    gc_min = values["GC_MIN_FRACTION"]
    gc_max = values["GC_MAX_FRACTION"]
    if not 0.0 <= gc_min < gc_max <= 1.0:
        raise SystemExit(
            "GC fraction must satisfy 0 <= GC_MIN_FRACTION < GC_MAX_FRACTION <= 1, "
            "got %s, %s" % (gc_min, gc_max)
        )
    if values["INSERT_SIZE"] <= 0:
        raise SystemExit("INSERT_SIZE must be a positive integer")


def render_value(name, value):
    if PARAM_TYPES[name] is float:
        return repr(float(value))
    return str(int(value))


def apply_replacements(text, updates):
    def repl(match):
        name = match.group("name")
        if name not in updates:
            return match.group(0)
        return "%s = %s" % (name, render_value(name, updates[name]))

    return ASSIGN_RE.sub(repl, text)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Patch six hard-coded constants in CNVkit params.py."
    )
    parser.add_argument("--file_path", type=Path, default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--force_rewrite",
        type=str2bool,
        default=None,
        help="False leaves the file unchanged. True is accepted and is not required.",
    )
    parser.add_argument(
        "--reference_auto_model",
        type=str2bool,
        default=False,
        help="Deprecated. Accepted for old commands. reference.py is not edited.",
    )
    for name, value_type in PARAM_TYPES.items():
        parser.add_argument("--%s" % name, type=value_type, default=None)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.reference_auto_model:
        print(
            "warning: --reference_auto_model is ignored. "
            "CNVkit 0.9.14 reconciles target and antitarget sex evidence in "
            "_reconcile_sex_guesses() and this image does not edit reference.py.",
            file=sys.stderr,
        )

    path = args.file_path if args.file_path is not None else default_params_path()
    if not path.is_file():
        raise SystemExit("params.py not found: %s" % path)

    original = path.read_text(encoding="utf-8")
    current = {
        name: parse_literal(name, raw) for name, raw in read_assignments(original).items()
    }
    requested = {
        name: getattr(args, name)
        for name in PARAM_TYPES
        if getattr(args, name) is not None
    }
    merged = dict(current)
    merged.update(requested)
    validate(merged)

    print("params file: %s" % path)
    for name in PARAM_TYPES:
        if name in requested and requested[name] != current[name]:
            print("%s: %s -> %s" % (name, current[name], requested[name]))
        else:
            print("%s: %s" % (name, current[name]))

    if not requested or args.force_rewrite is False:
        if args.force_rewrite is False and requested:
            print("force_rewrite is False; file not modified")
        return 0

    updated = apply_replacements(original, requested)
    if "REGISTERED_BUILDS" in original and "REGISTERED_BUILDS" not in updated:
        raise SystemExit("refusing to drop skgenome.genomebuild configuration")
    compile(updated, str(path), "exec")

    if updated == original:
        print("values already match; file not modified")
        return 0

    diff = "".join(
        difflib.unified_diff(
            original.splitlines(keepends=True),
            updated.splitlines(keepends=True),
            fromfile=str(path),
            tofile=str(path) + " (patched)",
        )
    )
    sys.stdout.write(diff)

    if args.dry_run:
        print("dry-run; file not modified")
        return 0

    backup = path.with_name("params.py.orig")
    if not backup.exists():
        shutil.copy2(path, backup)
    path.write_text(updated, encoding="utf-8")

    written = {
        name: parse_literal(name, raw)
        for name, raw in read_assignments(path.read_text(encoding="utf-8")).items()
    }
    for name, value in requested.items():
        if written[name] != value:
            raise SystemExit(
                "verification failed for %s: file has %s, expected %s"
                % (name, written[name], value)
            )
    print("patched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
