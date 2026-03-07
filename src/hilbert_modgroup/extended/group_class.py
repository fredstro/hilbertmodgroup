import logging

import sage
from sage.all import Integer
from sage.arith.misc import divisors
from sage.categories.groups import Groups
from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic
from sage.matrix.constructor import Matrix, matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modular.cusps_nf import NFCusps_ideal_reps_for_levelN, list_of_representatives
from sage.modular.modsym.p1list_nf import psi
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import QuadraticField

from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement__class

from .cusp import (
    NFCusp_wrt_lattice_ideal,
    fundamental_unit_generator,
    totally_positive_unit_group_generators,
)
from .group_element import ExtendedHilbertModularGroupElement

logger = logging.getLogger(__name__)


def is_ExtendedHilbertModularGroup(x) -> bool:
    """
    Return `True` if ``x`` is an instance of an ExtendedHilbertModularGroup

    INPUT:

    - ``x`` -- something to test if it is an Extended Hilbert modular group or not

    OUTPUT:

    - boolean

    EXAMPLES::

        sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
        sage: from hilbert_modgroup.extended.all import is_ExtendedHilbertModularGroup
        sage: is_ExtendedHilbertModularGroup(1)
        False
        sage: H = ExtendedHilbertModularGroup(5)
        sage: is_ExtendedHilbertModularGroup(H)
        True
        sage: K3.<a> = QuadraticField(3)
        sage: lattice_ideal = K3.fractional_ideal(5)
        sage: level_ideal = K3.fractional_ideal(7)
        sage: H = ExtendedHilbertModularGroup(K3, lattice_ideal, level_ideal)
        sage: is_ExtendedHilbertModularGroup(H)
        True
    """
    return isinstance(x, ExtendedHilbertModularGroup_class)


def ExtendedHilbertModularGroup(number_field, lattice_ideal=None, level_ideal=None, tp_units=True):
    r"""
    Create the Extended Hilbert modular group PGL_2^+(Fractional_ideal(1) \oplus lattice_ideal, level_ideal)
    (or PSL_2(Fractional_ideal(1) \oplus lattice_ideal, level_ideal)), of order 2, consisting of matrices of
    determinant in U_K^+ (or 1).

    INPUT:

    - ``number_field`` (NumberField)       -- a totally real field.
    - ``lattice_ideal`` (NumberFieldIdeal)  -- an ideal in totally real field.
    - ``level_ideal`` (NumberFieldIdeal)    -- an ideal in the same field.
    - ``tp_units `` (bool)  - ``True`` if you want the determinant in U_K^+ and ``False`` if 1.

    EXAMPLES::

        sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
        sage: number_field = 5
        sage: ExtendedHilbertModularGroup(number_field)
        Hilbert modular group PGL_2^+(Fractional ideal (1) + Fractional ideal (1), Fractional ideal (1)) of order
        2 over Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? consisting of
        matrices of determinant in U_K^+.
        sage: K1.<a> = QuadraticField(2)
        sage: lattice_ideal = K1.fractional_ideal(3)
        sage: ExtendedHilbertModularGroup(K1, lattice_ideal)
        Hilbert modular group PGL_2^+(Fractional ideal (1) + Fractional ideal (3), Fractional ideal (1)) of
        order 2 over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? consisting
        of matrices of determinant in U_K^+.
        sage: level_ideal = K1.fractional_ideal(5)
        sage: ExtendedHilbertModularGroup(K1, lattice_ideal, level_ideal)
        Hilbert modular group PGL_2^+(Fractional ideal (1) + Fractional ideal (3), Fractional ideal (5)) of
        order 2 over Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? consisting
        of matrices of determinant in U_K^+.
        sage: ExtendedHilbertModularGroup(K1, lattice_ideal, level_ideal, tp_units = False)
        Extended Hilbert modular group PSL_2(Fractional ideal (1) +  lattice_ideal, level_ideal) of order 2 over
        Number Field in a with defining polynomial x^2 - 2 with a = 1.414213562373095? consisting of matrices
        of determinant 1.

    """
    if isinstance(number_field, int | Integer):
        if number_field > 0:
            number_field = QuadraticField(number_field)
        else:
            raise ValueError("The input must be a positive integers")
    if (
        isinstance(number_field, sage.rings.number_field.number_field_base.NumberField)
        and not number_field.is_totally_real()
    ):
        raise ValueError("The input must be a totally real Number Field")
    if not lattice_ideal:
        lattice_ideal = number_field.fractional_ideal(1)
    if not level_ideal:
        level_ideal = number_field.fractional_ideal(1)
    Integer(2)
    if not (
        lattice_ideal.number_field() == number_field
        and level_ideal.number_field() == number_field
        and lattice_ideal.is_coprime(level_ideal)
        and level_ideal.is_coprime(number_field.different())
    ):
        raise ValueError(
            "The lattice ideal and level ideal must be an ideal of the same number field and "
            "coprime to each other and level ideal is coprime to different ideal of the number field"
        )
    if tp_units:
        name = (
            f"Hilbert modular group PGL_2^+({number_field.fractional_ideal(1)} + {lattice_ideal}, "
            f"{level_ideal}) of order 2 over {number_field} consisting of matrices of determinant in U_K^+."
        )
    else:
        name = (
            f"Extended Hilbert modular group PSL_2({number_field.fractional_ideal(1)} +  "
            f"lattice_ideal, level_ideal) of order 2 over {number_field} consisting of "
            f"matrices of determinant 1."
        )
    ltx = r"PGL_2^+[\mathcal{{O}}_K + \mathfrak{{a}}]"
    return ExtendedHilbertModularGroup_class(
        number_field,
        lattice_ideal=lattice_ideal,
        level_ideal=level_ideal,
        tp_units=tp_units,
        sage_name=name,
        latex_string=ltx,
    )


class ExtendedHilbertModularGroup_class(LinearMatrixGroup_generic):
    r"""
    Class for Extended Hilbert modular groups, here defined as either P[OK + I] (default) with determinants
    that are totally positive units in OK.
    """

    Element = ExtendedHilbertModularGroupElement

    def __init__(self, number_field, lattice_ideal, level_ideal, tp_units, sage_name, latex_string):
        r"""
             Init an Extended Hilbert modular group of the form PGL2^+[OK+lattice_ideal].

            INPUT:

            - ``lattice_ideal`` - NumberFieldFractionalIdeal
            - ``sage_name`` - string
            - ``latex_string`` - string

        EXAMPLES::

            sage: from hilbert_modgroup.extended.group_class import (
            ....:     ExtendedHilbertModularGroup_class
            ....: )
            sage: K5.<a> = QuadraticField(5)
            sage: number_field = K5
            sage: lattice_ideal = K5.different()
            sage: level_ideal = K5.fractional_ideal(7)
            sage: name = (
            ....:     f"Hilbert modular group PGL_2^+("
            ....:     f"{number_field.fractional_ideal(1)} + {lattice_ideal}, "
            ....:     f"{level_ideal}) of order 2 over {number_field} "
            ....:     f"consisting of matrices of determinant in U_K^+."
            ....: )
            sage: ltx = "PGL2^+[OK + lattice_ideal]"
            sage: G = ExtendedHilbertModularGroup_class(
            ....:     number_field = number_field,
            ....:     lattice_ideal = lattice_ideal,
            ....:     level_ideal = level_ideal,
            ....:     tp_units = True,
            ....:     sage_name = name,
            ....:     latex_string = ltx,
            ....: )
            sage: G
            Hilbert modular group PGL_2^+(Fractional ideal (1) + Fractional ideal (-a),
            Fractional ideal (7)) of order 2 over Number Field in a with defining
            polynomial x^2 - 5 with a = 2.236067977499790? consisting of matrices of
            determinant in U_K^+.
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: TestSuite(H1).run()
            sage: H1(1)
            [1 0]
            [0 1]
            sage: H1([1, 1, 0, 1])
            [1 1]
            [0 1]
            sage: H1([1, H1.OK().gens()[0], 0, 1])
            [          1 1/2*a + 1/2]
            [          0           1]

        """
        self._number_field = number_field
        self._lattice_ideal = lattice_ideal
        self._level_ideal = level_ideal
        self._OK = number_field.ring_of_integers()
        self._tp_units = tp_units
        # Instance data related to cusps
        self._ncusps = None
        self._cusps = []
        self._ideal_cusp_representatives = []
        self._cusp_normalizing_maps = {}
        self._cusp_normalizing_maps_inverse = {}
        super().__init__(
            degree=Integer(2),
            base_ring=number_field,
            special=False,
            sage_name=sage_name,
            latex_string=latex_string,
            category=Groups().Infinite(),
            invariant_form=None,
        )

    def number_field(self):
        """
        Return the number field associated to self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: H.number_field()
            Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
        """
        return self._number_field

    def lattice_ideal(self):
        """
        Return the lattice_ideal associated to self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: H = ExtendedHilbertModularGroup(5)
            sage: H.lattice_ideal()
            Fractional ideal (1)
        """
        return self._lattice_ideal

    def level_ideal(self):
        """
        Return the level_ideal associated to self.

        EXAMPLES::
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K = QuadraticField(3)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: H.level_ideal()
            Fractional ideal (1)
            sage: level_ideal = K.fractional_ideal(7)
            sage: H = ExtendedHilbertModularGroup(K, level_ideal = level_ideal)
            sage: H.level_ideal()
            Fractional ideal (7)

        """
        return self._level_ideal

    def tp_units(self):
        """
        Return the tp_units value associated to self.

        EXAMPLES::
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K = QuadraticField(3)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: H.tp_units()
            True
            sage: tp_units = False
            sage: H = ExtendedHilbertModularGroup(K, tp_units = False)
            sage: H.tp_units()
            False
        """
        return self._tp_units

    def OK(self):
        """
        Return the ring of integers associated with self.

        EXAMPLES::
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: H = ExtendedHilbertModularGroup(5)
            sage: H.OK()
            Maximal Order generated by 1/2*a + 1/2 in Number Field in a with defining polynomial x^2 - 5
            with a = 2.236067977499790?
        """
        return self._OK

    def __contains__(self, x):
        r"""
        Return whether ``x`` is an element of ``self``.

        This checks whether ``x`` can be converted into an element of this
        extended Hilbert modular group, i.e., whether it defines a valid
        2x2 matrix satisfying the group's determinant and ideal conditions.

        INPUT:

        - ``x`` -- any object that can be used to construct an
          :class:`ExtendedHilbertModularGroupElement`, e.g., a list of
          four elements or an existing group element

        OUTPUT:

        - boolean

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K5.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K5)
            sage: [1, 0, 0, 1] in H
            True
            sage: x, y = H.OK().gens()
            sage: ExtendedHilbertModularGroupElement(H, [1, x, 0, 1]) in H
            True
            sage: [1, 2, 3, 4] in H
            False
        """
        try:
            ExtendedHilbertModularGroupElement(self, x)
            return True
        except TypeError:
            return False

    def ideal(self, x: tuple):
        """
        Return the ideal generated by the tuple x given by x[0]O_K +x[1]lattice_ideal^{-1}.

        EXAMPLES::
           sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
           sage: K2.<a> = QuadraticField(2)
           sage: lattice_ideal = K2.fractional_ideal(3)
           sage: level_ideal = K2.fractional_ideal(1)
           sage: H2 = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
           sage: H2.ideal_cusp_representatives()
           [Fractional ideal (1)]
           sage: H2.ideal((2, 3))
           Fractional ideal (1)
           sage: H2.ideal((2, a))
           Fractional ideal (-1/3*a)
        """
        lattice_ideal = self.lattice_ideal()
        return NFCusp_wrt_lattice_ideal(lattice_ideal, x[0], x[1]).ideal()

    def generators(self):
        r"""
        Return a list of generators of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K1.<a> = QuadraticField(2)
            sage: H = ExtendedHilbertModularGroup(K1)
            sage: H.generators()
            [
            [1 1]  [1 a]  [1 0]  [1 0]  [2*a + 3       0]
            [0 1], [0 1], [1 1], [a 1], [      0       1]
            ]
            sage: lattice_ideal = K1.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K1, lattice_ideal)
            sage: H.generators()
            [
            [  1 1/3]  [    1 1/3*a]  [1 0]  [  1   0]  [2*a + 3       0]
            [  0   1], [    0     1], [3 1], [3*a   1], [      0       1]
            ]
            sage: lattice_ideal = K1.fractional_ideal(1)
            sage: level_ideal = K1.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K1, lattice_ideal, level_ideal)
            sage: H.generators()
            [
            [-a + 1      0]  [-1  0]  [a + 1     0]  [ a -1]  [-a -1]
            [     3 -a - 1], [ 3 -1], [    3 a - 1], [ 3 -a], [ 3  a],
            <BLANKLINE>
            [-a - 1      0]  [a - 1     0]  [1 1]  [1 a]  [1 0]  [  1   0]
            [     3 -a + 1], [    3 a + 1], [0 1], [0 1], [3 1], [3*a   1],
            <BLANKLINE>
            [2*a + 3       0]
            [      0       1]
            ]
            sage: H = ExtendedHilbertModularGroup(K1, lattice_ideal, level_ideal, tp_units = False)
            sage: H.generators()
            [
            [-a + 1      0]  [-1  0]  [a + 1     0]  [ a -1]  [-a -1]
            [     3 -a - 1], [ 3 -1], [    3 a - 1], [ 3 -a], [ 3  a],
            <BLANKLINE>
            [-a - 1      0]  [a - 1     0]  [1 1]  [1 a]  [1 0]  [  1   0]
            [     3 -a + 1], [    3 a + 1], [0 1], [0 1], [3 1], [3*a   1]
            ]

        """
        gens = []
        tp_units = self.tp_units()
        lattice_ideal = self.lattice_ideal()
        level_ideal = self.level_ideal()
        number_field = self.number_field()
        Lreps = list_of_representatives(level_ideal)
        for d in level_ideal.residues():
            if d != 0 and d != 1 and number_field.fractional_ideal(d).is_coprime(level_ideal):
                Lds = [
                    P * lattice_ideal * level_ideal
                    for P in Lreps
                    if (P * lattice_ideal * level_ideal).is_principal()
                ]
                C = Lds[0]
                c = (C).gens_reduced()[0]
                A1 = c * (lattice_ideal.inverse())
                A2 = number_field.fractional_ideal(d)
                r = A1.element_1_mod(A2)
                b = -r / c
                a = (1 - r) / d
                gens.append(self.create_element(a, b, c, d))
        for x in self.lattice_ideal().inverse().basis():
            gens.append(self.T(x))
        for x in (self.lattice_ideal() * self.level_ideal()).basis():
            gens.append(self.L(x))
        if tp_units:
            tpunit_gen = totally_positive_unit_group_generators(number_field)
            for x in tpunit_gen:
                gens.append(self.E(x))
        return gens

    @cached_method
    def S(self):
        """
        Return the element S = ( 0 & -1 // 1 & 0 ) of self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K1.<a> = QuadraticField(2)
            sage: H = ExtendedHilbertModularGroup(K1)
            sage: H.S()
            [ 0 -1]
            [ 1  0]
        """
        return self([0, -1, 1, 0])

    @cached_method
    def T(self, a=1):
        """
        Return the element T^a = ( 1 & a // 0 & 1 ) of self.

        INPUT:

        - ``a`` -- integer in number field (default=1)

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: H.T()
            [1 1]
            [0 1]
            sage: u0, u1 = K2.unit_group().gens()
            sage: H.T(u0)
            [ 1 -1]
            [ 0  1]
            sage: H.T(u0*u1)
            [          1 1/2*a - 1/2]
            [          0           1]
        """
        return self([1, a, 0, 1])

    @cached_method
    def L(self, a):
        """
        Return the element L=( 1 & 0 // a & 1 ) of self.

        INPUT:

        - ``a`` -- integer in number field

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: H = ExtendedHilbertModularGroup(5)
            sage: H.L(1)
            [1 0]
            [1 1]
            sage: u0,u1 = H.number_field().unit_group().gens()
            sage: H.L(u0)
            [ 1 0]
            [-1 1]
            sage: H.L(u0*u1)
            [          1           0]
            [1/2*a - 1/2           1]

        """
        return self([1, 0, a, 1])

    @cached_method
    def E(self, u):
        """
        Return the element U =( u & 0 // 0 & 1 ) of self.

        INPUT:
        When tp_unit is True then u must be an element of  totally positive group.
        Else u is an element of unit group.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: H5.E(1)
            [1 0]
            [0 1]
            sage: u0, u1 = H5.number_field().unit_group().gens()
            sage: H5.E(u0*u1)
            Traceback (most recent call last):
            ...
            ValueError: u must be totally positive unit.
            sage: H5 = ExtendedHilbertModularGroup(K5, tp_units = False)
            sage: H5.E(u0*u1)
            [1/2*a - 1/2           0]
            [          0 1/2*a + 1/2]

        """
        tp_units = self.tp_units()
        K = self.number_field()
        if tp_units:
            if not K(u).is_totally_positive():
                raise ValueError("u must be totally positive unit.")
            return self([u, 0, 0, 1])
        else:
            return self([u, 0, 0, u**-1])

    def create_element(self, a, b, c, d):
        r"""
        Return an element of Extended Hilbert Modular Group.

        INPUT:

        - a, b, c, d -- matrix entries.


        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: H = ExtendedHilbertModularGroup(5)
            sage: A = H.create_element(1, 1, 0, 1)
            sage: A in H
            True
        """
        return self([a, b, c, d])

    def random_element(self, matrix_type=None, **kwds):
        r"""
        Return a 'random' element of this Extended Hilbert Modular Group.

        INPUT:

        - ``mode`` -- one of {'Lower', 'Upper', 'unit'} or None (default)
        - ``kwds`` -- passed to the random element generators

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: H = ExtendedHilbertModularGroup(5)
            sage: A = H.random_element()
            sage: A in H
            True
        """
        x = kwds.pop("x", None)
        y = kwds.pop("y", None)
        a = self.lattice_ideal().inverse().random_element(**kwds)
        b = (self.lattice_ideal() * self.level_ideal()).random_element(**kwds)
        K = self.number_field()
        if x is None:
            x = -5
        if y is None:
            y = 5
        if self.tp_units():
            gens = totally_positive_unit_group_generators(K)
            exponents = [ZZ.random_element(x, y) for _ in gens]
            u = prod(g**e for g, e in zip(gens, exponents, strict=False))
        else:
            gens = fundamental_unit_generator(K)
            exponents = [ZZ.random_element(x, y) for _ in gens]
            u = prod(g**e for g, e in zip(gens, exponents, strict=False))
        if matrix_type == "Lower":
            return self(self.L(b))

        if matrix_type == "Upper":
            return self(self.T(a))

        if matrix_type == "unit":
            return self(self.E(u))

        return self(self.E(u) * self.T(a) * self.L(b))

    @cached_method
    def cusps(self):
        """
        Return the inequivalent cusps of self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K5.<a> = QuadraticField(5)
            sage: lattice_ideal = K5.fractional_ideal(1)
            sage: level_ideal = K5.fractional_ideal(3)
            sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal, level_ideal)
            sage: H5.cusps()
            [Cusp [0: 1] of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? with respect to  lattice_ideal,
             Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? with respect to lattice_ideal]
            sage: lattice_ideal = K5.fractional_ideal(7)
            sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal, level_ideal)
            sage: H5.cusps()
            [Cusp [0: 7] of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? with respect to  lattice_ideal,
             Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? with respect to lattice_ideal]
        """
        K = self.number_field()
        lattice_ideal = self.lattice_ideal()
        N = self.level_ideal()
        L = NFCusps_ideal_reps_for_levelN(N, nlists=3)
        Laux = L[1] + L[2]
        Lreps = self.ideal_cusp_representatives()
        Lcusps = []
        for A in L[0]:
            if A.is_trivial():
                B = K.ideal(1)
                g = 1
            else:
                Lbs = [P for P in Laux if (P * A).is_principal()]
                B = Lbs[0]
                g = (A * B).gens_reduced()[0]
            for d in divisors(N):
                Lds = [
                    P
                    for P in Laux
                    if (P * lattice_ideal * d * A).is_principal() and P.is_coprime(B)
                ]
                deltap = Lds[0]
                c = (deltap * lattice_ideal * d * A).gens_reduced()[0]
                I = d + N / d
                if c.is_one() and I.is_trivial():
                    Lcusps.append(NFCusp_wrt_lattice_ideal(lattice_ideal, 0, 1, lreps=Lreps))
                else:
                    if self.tp_units():
                        gens = totally_positive_unit_group_generators(K)
                    else:
                        u = K.unit_group().gens_values()
                        gens = [t**2 for t in u]
                    for b in I.invertible_residues_mod(gens):
                        M = d.prime_to_idealM_part(I)
                        deltAM = deltap * lattice_ideal * A * M
                        u = (B * deltAM).element_1_mod(I)
                        v = (I * B).element_1_mod(deltAM)
                        newb = u * b + v
                        A1 = c * (lattice_ideal.inverse()) * (A.inverse())
                        A2 = newb * B.inverse()
                        r = A2.element_1_mod(A1)
                        a1 = (r / newb) * g
                        -(1 - r) / c * g
                        Lcusps.append(NFCusp_wrt_lattice_ideal(lattice_ideal, a1, c, lreps=Lreps))
        cusp = NFCusp_wrt_lattice_ideal(self.lattice_ideal(), 1, 0)
        for c in Lcusps:
            t = cusp.is_Gamma0_wrt_lattice_ideal_equivalent(c, self.level_ideal())
            if t:
                idx = Lcusps.index(c)
                Lcusps[idx] = NFCusp_wrt_lattice_ideal(self.lattice_ideal(), 1, 0, lreps=Lreps)
                break
        self._cusps = Lcusps
        return self._cusps

    def ncusps(self):
        """
        Return the number of cusp associated to self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K2.<a> = QuadraticField(2)
            sage: lattice_ideal = K2.fractional_ideal(7)
            sage: level_ideal = K2.fractional_ideal(83)
            sage: H2 = ExtendedHilbertModularGroup(K2, lattice_ideal = lattice_ideal, level_ideal = level_ideal)
            sage: H2.ncusps()
            2
            sage: level_ideal = K2.fractional_ideal(9)
            sage: H2 = ExtendedHilbertModularGroup(K2, lattice_ideal = lattice_ideal, level_ideal = level_ideal)
            sage: H2.ncusps()
            4
        """
        K = self.number_field()
        N = self.level_ideal()
        if self.tp_units():
            gens = totally_positive_unit_group_generators(K)
        else:
            u = K.unit_group().gens_values()
            gens = [K(t**2) for t in u]
        s = sum([len((d + N / d).invertible_residues_mod(gens)) for d in divisors(N)])
        self._ncusps = s * K.class_number()
        return self._ncusps

    def ideal_cusp_representatives(self):
        r"""
        Return a list of ideals corresponding to cusp representatives, i.e.
        ideal representatives of ideal classes.

        Note: We choose an ideal of smallest norm in each class.
            If the ideal given by sage is already minimal we return this.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K2.<a> = QuadraticField(2)
            sage: lattice_ideal = K2.fractional_ideal(1)
            sage: level_ideal = K2.fractional_ideal(1)
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
            sage: H.ideal_cusp_representatives()
            [Fractional ideal (1)]
            sage: K10.<a> = QuadraticField(10)
            sage: lattice_ideal = K10.fractional_ideal(1)
            sage: level_ideal = K10.fractional_ideal(1)
            sage: H10 = ExtendedHilbertModularGroup(K10, lattice_ideal, level_ideal)
            sage: H10.ideal_cusp_representatives()
            [Fractional ideal (1), Fractional ideal (2, a)]
            sage: lattice_ideal = K10.fractional_ideal(2, a)
            sage: H10 = ExtendedHilbertModularGroup(K10, lattice_ideal, level_ideal)
            sage: H10.ideal_cusp_representatives()
            [Fractional ideal (1), Fractional ideal (2, a)]
            sage: level_ideal = K10.fractional_ideal(3)
            sage: H10 = ExtendedHilbertModularGroup(K10, lattice_ideal, level_ideal)
            sage: H10.ideal_cusp_representatives()
            [Fractional ideal (1), Fractional ideal (5, a)]


        """
        if not self._ideal_cusp_representatives:
            self._ideal_cusp_representatives = []

            def _find_equivalent_ideal_of_minimal_norm(c, ncpi):
                for a in self.number_field().ideals_of_bdd_norm(c.norm() - 1).items():
                    for ideala in a[1]:
                        if (ideala * c**-1).is_principal() and ideala.is_coprime(ncpi):
                            if c.norm() <= ideala.norm():
                                return c
                            return ideala
                return c

            lattice_ideal = self.lattice_ideal()
            N = self.level_ideal()
            if self.number_field().fractional_ideal(1) == N:
                ncpi = self.number_field().fractional_ideal(1)
            else:
                ncpi = lattice_ideal * N
            Lreps = list_of_representatives(ncpi)
            for c in Lreps:
                c = _find_equivalent_ideal_of_minimal_norm(c, ncpi)
                self._ideal_cusp_representatives.append(c)
            self._ideal_cusp_representatives.sort(key=lambda x: x.norm())
        return self._ideal_cusp_representatives

    def cusp_representative(self, cusp, return_map=False):
        r"""
        Return a representative cusp and optionally a corresponding map.

        INPUT:
        - ``cusp`` -- cusp
        - ``return_map`` -- bool (default: False)
            Set to True to also return the map giving the equivalence.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K10.<a> = QuadraticField(10)
            sage: lattice_ideal = K10.different()
            sage: level_ideal = K10.fractional_ideal(7)
            sage: H10 = ExtendedHilbertModularGroup(K10, lattice_ideal, level_ideal)
            sage: from hilbert_modgroup.extended.all import NFCusp_wrt_lattice_ideal
            sage: c = NFCusp_wrt_lattice_ideal(lattice_ideal, 0, 1)
            sage: H10.cusp_representative(c)
            Cusp [0: -2*a] of Number Field in a with defining polynomial x^2 - 10 with
            a = 3.162277660168380? with respect to  lattice_ideal
            sage: a = H10.number_field().gen()
            sage: x, y = 3*a-10, a-4
            sage: c = NFCusp_wrt_lattice_ideal(lattice_ideal, x, y)
            sage: H10.cusp_representative(c)
            Cusp [-2*a - 1: -14*a + 80] of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380? with respect to  lattice_ideal
            sage: K5.<a> = QuadraticField(5)
            sage: lattice_ideal = K5.fractional_ideal(2)
            sage: level_ideal = K5.fractional_ideal(3)
            sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal, level_ideal, False)
            sage: c = NFCusp_wrt_lattice_ideal(lattice_ideal, 2, a)
            sage: H5.cusp_representative(c)
            Cusp [0: 2] of Number Field in a with defining polynomial x^2 - 5 with
            a = 2.236067977499790? with respect to  lattice_ideal

        """
        tp_units = self.tp_units()
        for c in self.cusps():
            if return_map:
                t, B = cusp.is_Gamma0_wrt_lattice_ideal_equivalent(
                    c, self.level_ideal(), tp_units=tp_units, Transformation=True
                )
                if t:
                    return c, self(B)
            elif cusp.is_Gamma0_wrt_lattice_ideal_equivalent(
                c, self.level_ideal(), tp_units=tp_units
            ):
                return c
        raise ArithmeticError(f"Could not find cusp representative for {cusp}")

    def cusp_normalizing_map(self, cusp, inverse=False, check=False):
        r"""
        Given a cusp (a:c) Return a matrix A = [[ a ,b ], [c , d]] in SL(2,K) such that
        A(Infinity)=(a:c).

        INPUT:
        - ``cusp`` -- Instance of NFCusp
        - ``inverse`` -- bool (default: False) set to True to return the inverse map
        - ``check`` -- bool (default: False) set to True to check the result

        If inverse = True then return A^-1

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K3.<a> = QuadraticField(3)
            sage: lattice_ideal = K3.different()
            sage: level_ideal = K3.fractional_ideal(1)
            sage: H3 = ExtendedHilbertModularGroup(K3, lattice_ideal, level_ideal)
            sage: H3.cusp_normalizing_map(H3.cusps()[0])
            [1 0]
            [0 1]
            sage: K5.<a> = QuadraticField(5)
            sage: lattice_ideal = K5.fractional_ideal(1)
            sage: level_ideal = K5.fractional_ideal(3)
            sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal, level_ideal)
            sage: H5.ncusps()
            2
            sage: H5.cusps()[1]
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with
            a = 2.236067977499790? with respect to lattice_ideal
            sage: H5.cusps()[0]
            Cusp [0: 1] of Number Field in a with defining polynomial x^2 - 5 with
            a = 2.236067977499790? with respect to  lattice_ideal
            sage: H5.cusp_normalizing_map(H5.cusps()[0])
             [ 0 -1]
             [ 1  0]
        """
        base_nf = self.number_field()
        if not isinstance(cusp, NFCusp_wrt_lattice_ideal) or cusp.number_field() != base_nf:
            raise ValueError(f"Input should be a NF cusp defined over {base_nf}!")
        ca, cb = (cusp.numerator(), cusp.denominator())
        if (ca, cb) not in self._cusp_normalizing_maps:
            a, b, c, d = cusp.GHmatrix_wrt_lattice_ideal()
            det = a * d - b * c
            A = Matrix(self.number_field(), 2, 2, [a, b / det, c, d / det])
            if check:
                infinity_cusp = NFCusp_wrt_lattice_ideal(self.lattice_ideal(), 1, 0)
                if infinity_cusp.apply(A.list()) != cusp or A.det() != 1:
                    msg = f"Did not get correct normalizing map A={A} to cusp: {cusp}"
                    raise ArithmeticError(msg)
            logger.debug("A=%s", A)
            logger.debug("A.det()=%s", A.det().complex_embeddings())
            self._cusp_normalizing_maps_inverse[(ca, cb)] = A.inverse()
            self._cusp_normalizing_maps[(ca, cb)] = A
        if inverse:
            return self._cusp_normalizing_maps_inverse[(ca, cb)]
        else:
            return self._cusp_normalizing_maps[(ca, cb)]

    def apply_cusp_normalizing_map(self, cusp, z, inverse=False):
        """
        Apply the cusp normalizing map associated with the cusp to an element z.

        INPUT:
        - `cusp` - an instance of NFcusp_wrt_lattice_ideal
        - `z` - an element in
             - the base number field
             - the set for cusps wrt lattice_ideal
             - in ComplexPlaneProductElement__class
        - `inverse` - set to True if applying the inverse map

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import NFCusp_wrt_lattice_ideal
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K10.<a> = QuadraticField(10)
            sage: lattice_ideal = K10.fractional_ideal(1)
            sage: level_ideal = K10.fractional_ideal(1)
            sage: H10 = ExtendedHilbertModularGroup(K10, lattice_ideal, level_ideal)
            sage: H10.apply_cusp_normalizing_map(H10.cusps()[1], 1.0)
            -1.30229210112909
            sage: z = NFCusp_wrt_lattice_ideal(lattice_ideal, 1)
            sage: H10.apply_cusp_normalizing_map(H10.cusps()[1], z)
            Cusp [a - 6: -5*a + 18] of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380? with respect to  lattice_ideal
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z = ComplexPlaneProductElement([CC(1,0),CC(1,0)]); z
            [1.00000000000000, 1.00000000000000]
            sage: H10.apply_cusp_normalizing_map(H10.cusps()[1], z)
            [-0.271492527522411, -1.30229210112909]
            sage: H10.apply_cusp_normalizing_map(H10.cusps()[1], 1).complex_embeddings()
            [-0.271492527522411, -1.30229210112909]
        """
        a, b, c, d = self.cusp_normalizing_map(cusp, inverse=inverse).list()
        if z == infinity:
            return a / c
        number_field = self.number_field()
        if isinstance(z, NFCusp_wrt_lattice_ideal) and z.number_field() == number_field:
            return z.apply([a, b, c, d])
        if z in number_field:
            return (a * z + b) / (c * z + d)
        if (
            isinstance(z, ComplexPlaneProductElement__class)
            and z.degree() == number_field.absolute_degree()
        ):
            return z.apply(matrix(2, 2, [a, b, c, d]))
        raise ValueError(f"Unsupported type for action with cusp normalizer! (z={z})")

    @cached_method
    def coset_matrices(self):
        """
        Return the list of right cosets representatives of GL_2^+(OK oplus ida, idn) in GL_2(OK oplus ida).

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K5.<a> = QuadraticField(5)
            sage: lattice_ideal = K5.fractional_ideal(1)
            sage: level_ideal = K5.fractional_ideal(3)
            sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal, level_ideal)
            sage: H5.coset_matrices()
            [
            [           0           -1]  [-1/2*a + 1/2            0]
            [           1 -3/2*a - 1/2], [           1 -1/2*a - 1/2],
            <BLANKLINE>
            [1/2*a + 1/2           0]  [ 0 -1]  [ 1 -1]  [ 0 -1]
            [          1 1/2*a - 1/2], [ 1 -a], [ 1  0], [ 1  a],
            <BLANKLINE>
            [-1/2*a - 1/2            0]  [1/2*a - 1/2           0]
            [           1 -1/2*a + 1/2], [          1 1/2*a + 1/2],
            <BLANKLINE>
            [          0          -1]  [1 0]
            [          1 3/2*a + 1/2], [0 1]
            ]
            sage: from sage.modular.modsym.p1list_nf import psi
            sage: psi(level_ideal)
            10
            sage: len(H5.coset_matrices())
            10
            sage: level_ideal = K5.fractional_ideal(7)
            sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal, level_ideal)
            sage: len(H5.coset_matrices()) == psi(level_ideal)
            True
        """
        N = self.level_ideal()
        K = self.number_field()
        lattice_ideal = self.lattice_ideal()
        H = ExtendedHilbertModularGroup(K, lattice_ideal)
        L = []
        for D in divisors(N):
            if (D * lattice_ideal).is_principal():
                Dp = K.fractional_ideal(1)
                c = (D * lattice_ideal * Dp).gens_reduced()[0]
            else:
                it = K.primes_of_degree_one_iter()
                Dp = next(it)
                while not Dp.is_coprime(N) or not (Dp * D * lattice_ideal).is_principal():
                    Dp = next(it)
                c = (D * lattice_ideal * Dp).gens_reduced()[0]
            I = D + N / D
            for r in (N / D).residues():
                if I.is_coprime(r):
                    M = D.prime_to_idealM_part(N / D)
                    u = (Dp * M).element_1_mod(N / D)
                    d = u * r + (1 - u)
                    if d.is_zero():
                        L.append(H.create_element(1, -1 / c, c, d))
                    else:
                        B = K.fractional_ideal(c * lattice_ideal.inverse()).element_1_mod(
                            K.fractional_ideal(d)
                        )
                        b = -B / c
                        a = (1 - B) / d
                        L.append(H.create_element(a, b, c, d))
        for x in L:
            if x in self:
                idx = L.index(x)
                L[idx] = H.create_element(1, 0, 0, 1)
                break
        if not len(L) == psi(N):
            raise ValueError("Condition is not satisfying. Check again")
        return L
