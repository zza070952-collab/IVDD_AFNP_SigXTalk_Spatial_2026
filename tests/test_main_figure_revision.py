"""Discover the revision module's technical tests in the existing CI job."""
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]/'scripts/authoring/main_figures_20261025'


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for folder in ['references', 'figures2_6']:
        location = ROOT/folder
        if location.is_dir():
            suite.addTests(unittest.TestLoader().discover(
                str(location), pattern='test_*.py', top_level_dir=str(location)))
    return suite
