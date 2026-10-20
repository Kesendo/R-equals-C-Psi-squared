"""Shared pytest configuration for the framework tests: the `slow` marker.

Tests marked slow run exact computations of minutes; they are skipped unless pytest is given --runslow."""
import pytest


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="run the tests marked slow")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: an exact run of minutes; skipped unless --runslow is given")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip = pytest.mark.skip(reason="slow; run with --runslow")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)
