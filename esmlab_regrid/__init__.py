#!/usr/bin/env python
"""Top-level module for esmlab-regrid."""

from pkg_resources import DistributionNotFound, get_distribution

from . import config
from .core import Regridder as regridder

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
