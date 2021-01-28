"""
Test the module
"""
from types import ModuleType
import spdelib


def test_src():
    assert isinstance(spdelib, ModuleType)
