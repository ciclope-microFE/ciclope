"""Computed Tomography to Finite Elements"""
__version__ = "2.0.1"

from .core import voxelFE
try:
    from .core import tetraFE
except ImportError:
    import warnings
    warnings.warn("tetraFE not available.", RuntimeWarning)

# from .utils import preprocess, postprocess
