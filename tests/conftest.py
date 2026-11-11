"""Shared execution of examples used by documentation and application checks."""

import runpy
from functools import cache

import pytest


@pytest.fixture(scope="session")
def run_example():
    """Execute each script once and retain its namespace for independent checks."""

    @cache
    def execute(path):
        return runpy.run_path(str(path), run_name="__main__")

    return execute


def pytest_addoption(parser):
    parser.addoption(
        "--without-real", action="store_true", help="collect core tests without importing semialg"
    )


def pytest_ignore_collect(collection_path, config):
    if (
        config.getoption("--without-real")
        and collection_path.name == "real"
        and collection_path.is_dir()
    ):
        return True
    return None


def pytest_collection_modifyitems(config, items):
    if config.getoption("--without-real"):
        skip = pytest.mark.skip(reason="requires the real extra")
        for item in items:
            if item.get_closest_marker("real"):
                item.add_marker(skip)
