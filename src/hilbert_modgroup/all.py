try:
    # Need to import some symbols explicitly for command line imports when using Passagemath
    from sage.all__sagemath_symbolics import Category
except ImportError:
    pass
from .hilbert_modular_group_class import HilbertModularGroup, is_HilbertModularGroup
from .pullback import HilbertPullback
from .upper_half_plane import UpperHalfPlaneProductElement, ComplexPlaneProductElement, \
    ComplexPlaneProduct, UpperHalfPlaneProduct
