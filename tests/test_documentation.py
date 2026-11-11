"""Execute documented workflows and validate links into the source tree."""

import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest

REAL_EXAMPLES = {
    row["file"]
    for row in json.loads(
        (Path(__file__).resolve().parents[1] / "examples/catalog.json").read_text(encoding="utf-8")
    )
    if row.get("profile") == "real"
}
REAL_DOCUMENTS = {"docs/applications.md"}

ROOT = Path(__file__).resolve().parents[1]
DOCUMENTS = (ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md")), ROOT / "examples/README.md")


@pytest.mark.parametrize(
    "example",
    [
        pytest.param(p, marks=pytest.mark.real) if p.name in REAL_EXAMPLES else p
        for p in sorted((ROOT / "examples").glob("*.py"))
    ],
    ids=lambda p: p.stem,
)
def test_examples(example, run_example):
    run_example(example)


@pytest.mark.parametrize(
    "document",
    [
        pytest.param(p, marks=pytest.mark.real) if str(p.relative_to(ROOT)) in REAL_DOCUMENTS else p
        for p in DOCUMENTS
        if "```python" in p.read_text(encoding="utf-8")
    ],
    ids=lambda p: str(p.relative_to(ROOT)),
)
def test_python_blocks(document):
    text = document.read_text(encoding="utf-8")
    namespace = {}
    for match in re.finditer(r"^```python[^\S\n]*\n(.*?)^```[^\S\n]*$", text, re.M | re.S):
        # Preserve document line numbers when a snippet raises an exception.
        source = "\n" * text.count("\n", 0, match.start(1)) + match.group(1)
        exec(compile(source, str(document), "exec"), namespace)


@pytest.mark.parametrize("document", DOCUMENTS, ids=lambda p: str(p.relative_to(ROOT)))
def test_source_links(document):
    text = re.sub(r"```.*?```", "", document.read_text(encoding="utf-8"), flags=re.S)
    github = "https://github.com/BhuvaneshBhatt/alggeom/"
    for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)", text):
        if target.startswith(github):
            location = target[len(github) :].split("/", 2)
            if len(location) == 3 and location[:2] in (["blob", "main"], ["tree", "main"]):
                path = ROOT / unquote(urlsplit(location[2]).path)
                assert path.is_file() if location[0] == "blob" else path.is_dir(), target
        elif not urlsplit(target).scheme and not target.startswith("#"):
            path = document.parent / unquote(urlsplit(target).path)
            assert path.exists(), target
