#!/usr/bin/env python
"""Top-level module for esmlab-regrid."""


from . import config
from ._version import get_versions
from .core import Regridder as regridder

__version__ = get_versions()['version']
del get_versions
