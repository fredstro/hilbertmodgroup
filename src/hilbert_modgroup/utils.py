"""
General utility functions

"""
from sage.rings.real_mpfi import RealIntervalField


def upper(x, prec):
    r"""
    Return an upper bound to x given by precision prec.

    EXAMPLES::

        sage: from hilbert_modgroup.utils import upper
        sage: upper(0.1,8)
        0.11
        sage: upper(0.1,53)
        0.100000000000001
        sage: upper(RR.pi(),16)
        3.142
    """
    return RealIntervalField(prec)(x).upper()


def lower(x, prec):
    r"""
    Return a lower bound to x given by precision prec.

    EXAMPLES::

        sage: from hilbert_modgroup.utils import lower
        sage: lower(0.1,8)
        0.099
        sage: lower(0.1,30)
        0.099999999
        sage: lower(RR.pi(),16)
        3.141
    """
    return RealIntervalField(prec)(x).lower()
