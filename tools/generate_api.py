"""Generate the Markdown API reference from exported objects and their docstrings."""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import inspect
import pkgutil
from functools import cached_property
from pathlib import Path
from types import ModuleType

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "docs/api.md"
GITHUB = "https://github.com/BhuvaneshBhatt/alggeom/blob/main/"


def signature(value, name):
    """Render callable arguments without implementation-specific type representations."""
    if (
        inspect.isclass(value)
        and issubclass(value, BaseException)
        and "__init__" not in value.__dict__
    ):
        return f"class {name}({value.__bases__[0].__name__})"
    parameters = inspect.signature(value)
    parameters = parameters.replace(
        parameters=[
            p.replace(annotation=inspect.Parameter.empty) for p in parameters.parameters.values()
        ],
        return_annotation=inspect.Signature.empty,
    )
    prefix = "class " if inspect.isclass(value) else ""
    return f"{prefix}{name}{parameters}"


def documentation(value):
    """Convert docstring paragraphs to Markdown without arbitrary prose line breaks."""
    docstring = inspect.getdoc(value)
    if not docstring:
        raise ValueError(f"missing public docstring: {value}")
    return "\n\n".join(" ".join(p.splitlines()) for p in docstring.split("\n\n")).replace("``", "`")


def source_link(value):
    owner = value.__module__.split(".", 1)[0]
    package = importlib.import_module(owner)
    source = Path(package.__file__).resolve().parent.parent
    path = Path(inspect.getsourcefile(value)).resolve().relative_to(source)
    line = inspect.getsourcelines(value)[1]
    github = (
        GITHUB if owner == "alggeom" else f"https://github.com/BhuvaneshBhatt/{owner}/blob/main/"
    )
    return f"[Source]({github}src/{path.as_posix()}#L{line})"


def render_reference():
    """Document every public export, including certificate fields and inherited public members."""
    package = importlib.import_module("alggeom")
    modules = []
    modules.extend(
        importlib.import_module(info.name)
        for info in sorted(
            pkgutil.walk_packages(package.__path__, package.__name__ + "."),
            key=lambda info: info.name,
        )
        if not any(part.startswith("_") for part in info.name.split("."))
    )
    lines = [
        "# API reference",
        "",
        "This reference is generated from public signatures and docstrings by `python tools/generate_api.py`. Use `python tools/generate_api.py --check` to detect stale documentation.",
        "",
        "Start with the [geometry tutorial](geometry.md) for algebraic conventions, then [schemes and fibers](schemes.md), [components](components.md), [resolutions](resolutions.md), [curves](curves.md), [module arithmetic](modules.md), [local geometry](local_geometry.md), [finite fields](finite_fields.md) and [normalization](normalization.md), or the [decomposition guide](guide.md) for functional composition. See [result semantics](results.md) for decomposition certificates and the [coefficient-field tutorial](coefficient_fields.md) for implicit reconstruction. The decomposition `parameters` and `real` modules require the `real` extra. Generating the complete reference also requires that extra.",
        "",
    ]
    for module in modules:
        lines.extend([f"## `{module.__name__}`", "", documentation(module), ""])
        for name in module.__all__:
            value = getattr(module, name)
            if isinstance(value, ModuleType):
                continue
            lines.extend(
                [
                    f"### `{name}`",
                    "",
                    "```text",
                    signature(value, name),
                    "```",
                    "",
                    documentation(value),
                    "",
                    source_link(value),
                    "",
                ]
            )
            if inspect.isclass(value):
                if dataclasses.is_dataclass(value):
                    fields = ", ".join(f"`{field.name}`" for field in dataclasses.fields(value))
                    lines.extend([f"Stored fields: {fields}.", ""])
                members = {}
                for base in reversed(value.__mro__[:-1]):
                    if base.__module__.startswith("alggeom."):
                        members.update(base.__dict__)
                for member_name, member in members.items():
                    if member_name.startswith("_"):
                        continue
                    if isinstance(member, (property, cached_property)):
                        lines.extend(
                            [f"#### `{name}.{member_name}`", "", documentation(member), ""]
                        )
                    elif inspect.isfunction(member):
                        lines.extend(
                            [
                                f"#### `{name}.{member_name}`",
                                "",
                                "```text",
                                signature(member, member_name),
                                "```",
                                "",
                                documentation(member),
                                "",
                            ]
                        )
    return "\n".join(lines)


def main():
    """Write the generated reference or fail when its checked-in copy is stale."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = render_reference()
    if args.check:
        if REFERENCE.read_text(encoding="utf-8") != rendered:
            parser.exit(1, "docs/api.md is stale; run python tools/generate_api.py\n")
        print("Public API reference is current")
    else:
        REFERENCE.write_text(rendered, encoding="utf-8", newline="\n")
        print(f"Generated {REFERENCE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
