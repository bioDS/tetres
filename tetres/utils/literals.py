"""
Module contains Literals that are used for function input across the package.
"""
from typing import Literal

_MDS_TYPES = Literal["tsne"]
# TODO the dist is more a core thing that should be within a different module...
_DIST = Literal["rnni", "rf"]
