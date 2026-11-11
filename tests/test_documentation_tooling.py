"""Exercise released dependency selection and non-ASCII installed tutorials."""

import importlib.util
import sys
import tomllib
from pathlib import Path

import pytest
from packaging.requirements import Requirement

ROOT = Path(__file__).resolve().parents[1]


def _tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("profile", ["core", "real"])
@pytest.mark.parametrize("published", [False, True])
def test_backend_selection(monkeypatch, tmp_path, profile, published):
    installer = _tool("install_backends")
    for name in ("algroots", "semialg"):
        directory = tmp_path / name
        directory.mkdir()
        (directory / "pyproject.toml").write_text("[project]\n", encoding="utf-8")
    arguments = ["install", "--profile", profile]
    if published:
        monkeypatch.setattr(installer, "ROOT", tmp_path / "alggeom")
        (installer.ROOT / "ci").mkdir(parents=True)
        (installer.ROOT / "ci/backends.toml").write_bytes((ROOT / "ci/backends.toml").read_bytes())
        arguments.append("--published")
    else:
        arguments.extend(["--source-root", str(tmp_path)])
    monkeypatch.setattr(sys, "argv", arguments)
    calls = []
    monkeypatch.setattr(
        installer.subprocess, "run", lambda command, **kwargs: calls.append(command)
    )
    installer.main()
    names = ["algroots"] if profile == "core" else ["algroots", "semialg"]
    requirements = tomllib.loads((ROOT / "ci/backends.toml").read_text(encoding="utf-8"))
    expected = [
        requirements[name]["requirement"] if published else str((tmp_path / name).resolve())
        for name in names
    ]
    assert calls == [[sys.executable, "-m", "pip", "install", *expected]]


def test_backend_metadata():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    inputs = tomllib.loads((ROOT / "ci/backends.toml").read_text(encoding="utf-8"))
    core = {Requirement(value).name: Requirement(value) for value in project["dependencies"]}
    real = {
        Requirement(value).name: Requirement(value)
        for value in project["optional-dependencies"]["real"]
    }
    assert "semialg" not in core
    assert real["semialg"] == Requirement(inputs["semialg"]["requirement"])
    assert core["algroots"] == Requirement(inputs["algroots"]["requirement"])
    assert project["license"] == "GPL-3.0-only"


def test_utf8_tutorial(monkeypatch, tmp_path, capsys):
    checker = _tool("check_examples")
    (tmp_path / "examples").mkdir()
    (tmp_path / "docs").mkdir()
    (tmp_path / "examples/catalog.json").write_text("[]", encoding="utf-8")
    fence = chr(96) * 3
    (tmp_path / "README.md").write_text(
        f"# Rationals ℚ\n\n{fence}python\nassert 'ℚ' == chr(0x211A)\nprint('ℚ')\n{fence}\n",
        encoding="utf-8",
    )
    read_text = Path.read_text

    def windows_text(path, encoding=None, errors=None):
        return read_text(path, encoding=encoding or "cp1252", errors=errors)

    monkeypatch.setattr(Path, "read_text", windows_text)
    monkeypatch.setattr(checker, "ROOT", tmp_path)
    monkeypatch.setattr(checker, "installed", lambda: None)
    monkeypatch.setattr(sys, "argv", ["check", "--profile", "core"])
    monkeypatch.setenv("PYTHONUTF8", "0")
    checker.main()
    assert (
        "0 installed examples and 1 executable documentation blocks passed"
        in capsys.readouterr().out
    )
