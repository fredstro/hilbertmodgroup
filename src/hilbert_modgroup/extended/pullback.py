import logging

from sage.all import ZZ
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import floor
from sage.matrix.all import Matrix
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RealField
from sage.structure.sage_object import SageObject

from hilbert_modgroup.upper_half_plane import (
    ComplexPlaneProductElement__class,
    UpperHalfPlaneProductElement,
    UpperHalfPlaneProductElement__class,
)
from hilbert_modgroup.utils import lower, upper

from .cusp import (
    NFCusp_wrt_lattice_ideal,
    fundamental_unit_generator,
    totally_positive_unit_group_generators,
)
from .group_class import ExtendedHilbertModularGroup, ExtendedHilbertModularGroup_class
from .pullback_cython import (
    coordinates_to_ideal_elements,
    distance_to_cusp_eg,
    find_candidate_cusps,
    find_closest_cusp,
    lattice_elements_in_box,
)

log = logging.getLogger(__name__)


def stable_floor(x, eps=1e-12):
    r"""
    Deterministic replacement for ``floor(x + 1/2)``.

    Breaks ties consistently when ``x`` is near half-integers by
    subtracting a small epsilon before taking the floor.

    INPUT:

    - ``x`` -- a real number
    - ``eps`` -- a small positive real number (default: ``1e-12``)

    OUTPUT:

    An integer, equal to ``floor(x + 0.5 - eps)``.

    EXAMPLES::

        sage: from hilbert_modgroup.extended.pullback import stable_floor
        sage: stable_floor(1.3)
        1
        sage: stable_floor(1.5)
        1
        sage: stable_floor(1.5 + 1e-11)
        2
        sage: stable_floor(-0.7)
        -1
    """
    return floor(x + 0.5 - eps)


class ExtendedHilbertPullback(SageObject):
    def __init__(self, G):
        r"""
        Initialize the pullback for the extended Hilbert modular group ``G``.

        INPUT:

        - ``G`` -- an extended Hilbert modular group

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
            sage: K.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: P = ExtendedHilbertPullback(H)
            sage: P
            Pullback class for Hilbert modular group PGL_2^+(...) ... x^2 - 5 with a = 2.236067977499790?...

        """
        if not isinstance(G, ExtendedHilbertModularGroup_class):
            raise ValueError("Need a Extended Hilbert modular group")
        self._group = G

    def __eq__(self, other):
        r"""
        Check if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
            sage: K5.<a> = QuadraticField(5)
            sage: H1 = ExtendedHilbertModularGroup(K5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: P12 = ExtendedHilbertPullback(H1)
            sage: P1 == P12
            True
            sage: K3.<a> = QuadraticField(3)
            sage: H2 = ExtendedHilbertModularGroup(K3)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: P1 == P2
            False

        """
        return not (not isinstance(other, type(self)) or self._group != other._group)

    def __str__(self):
        r"""
        Return string representation of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
            sage: K.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: P = ExtendedHilbertPullback(H)
            sage: str(P)
            'Pullback class for Hilbert modular group PGL_2^+(...) ... x^2 - 5 with a = 2.236067977499790?...'

        """
        return f"Pullback class for {self._group}"

    def __repr__(self):
        r"""
        Return string representation of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
            sage: K.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: P = ExtendedHilbertPullback(H)
            sage: P
            Pullback class for Hilbert modular group PGL_2^+(...) ... x^2 - 5 with a = 2.236067977499790?...

        """
        return str(self)

    def group(self):
        r"""
        Return the group of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
            sage: K.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K)
            sage: P = ExtendedHilbertPullback(H)
            sage: P.group()
            Hilbert modular group PGL_2^+(...) ... x^2 - 5 with a = 2.236067977499790?...

        """
        return self._group

    def _check_upper_half_plane_element(self, z):
        r"""
        Check if z is an element of type UpperHalfPlaneProductElement__class
         of the correct degree.

        INPUT:
        - `z` - potential element of a product of upper half planes.

        EXAMPLES::
            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: P = ExtendedHilbertPullback(H)
            sage: P._check_upper_half_plane_element([1,1])
            Traceback (most recent call last):
            ...
            ValueError: Need ...UpperHalfPlaneProductElement__class of degree 2
            sage: z = UpperHalfPlaneProductElement([1+I,1+I])
            sage: P._check_upper_half_plane_element(z)
            True
            sage: z = UpperHalfPlaneProductElement([1+I,1+I,1+I])
            sage: P._check_upper_half_plane_element(z)
            Traceback (most recent call last):
            ...
            ValueError: Need ...UpperHalfPlaneProductElement__class of degree 2
            sage: z = ComplexPlaneProductElement([1+I,1+I])
            sage: P._check_upper_half_plane_element(z)
            Traceback (most recent call last):
            ...
            ValueError: Need ...UpperHalfPlaneProductElement__class of degree 2


        """
        if (
            not isinstance(z, UpperHalfPlaneProductElement__class)
            or z.degree() != self.group().number_field().degree()
        ):
            msg = (
                f"Need an element of type: "
                f"UpperHalfPlaneProductElement__class of degree "
                f"{self.group().number_field().degree()}"
            )
            raise ValueError(msg)
        return True

    @cached_method()
    def fundamental_units(self):
        r"""
        Return fundamental the units for the group of self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: P5.fundamental_units()
            [-1/2*a + 3/2]
            sage: K4.<a> = NumberField(x**3 - 36*x -1)
            sage: H4 = ExtendedHilbertModularGroup(K4, tp_units = True)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.fundamental_units()
            [a^2, a + 6]
            sage: H4 = ExtendedHilbertModularGroup(K4, tp_units = False)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.fundamental_units()
            [a, a + 6]
        """
        K = self.group().number_field()
        tp_units = self.group().tp_units()
        if tp_units:
            return totally_positive_unit_group_generators(K)
        else:
            return fundamental_unit_generator(K)

    @cached_method
    def basis_matrix_logarithmic_unit_lattice(self, prec=53):
        """
        Return the Basis matrix for the logarithmic unit lattice / +-1

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: P5.basis_matrix_logarithmic_unit_lattice()
            [ 0.962423650119207]
            [-0.962423650119207]
            sage: H5 = ExtendedHilbertModularGroup(K5, tp_units = False)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: P5.basis_matrix_logarithmic_unit_lattice()
            [ 0.481211825059603]
            [-0.481211825059603]
            sage: K3.<a> = QuadraticField(3)
            sage: H3 = ExtendedHilbertModularGroup(K3)
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: P3.basis_matrix_logarithmic_unit_lattice()
            [-1.31695789692482]
            [ 1.31695789692482]
            sage: H3 = ExtendedHilbertModularGroup(K3, tp_units = False)
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: P3.basis_matrix_logarithmic_unit_lattice()
            [-1.31695789692482]
            [ 1.31695789692482]

        """
        n = self.group().number_field().degree()
        entries = [
            [x.abs().log() for x in u.complex_embeddings(prec)] for u in self.fundamental_units()
        ]
        return matrix(RealField(prec), n - 1, n, entries).transpose()

    def Y(self, z, return_error_estimate=False):
        r"""
        Compute the coordinate of y=Im(z)/N(y)^(1/n)
        with respect to the logarithmic unit lattice.

        INPUT:

        - ``z`` -- element of type UpperHalfPlaneProductElement
        - ``return_error_estimate`` -- boolean (default=False) set to True
                                        to return a tuple including an estimate
                                        of the error.

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: u0, u1 = H5.number_field().unit_group().gens()
            sage: z = UpperHalfPlaneProductElement([CC(1, 1), CC(1, 1)])
            sage: P5.Y(z)
            (0.000000000000000)
            sage: P5.Y(u1*z)
            Traceback (most recent call last):
            ...
            ValueError: Need an element of type: UpperHalfPlaneProductElement__class of degree 2
            sage: P5.Y(u1**2*z)
            (1.00000000000000)
            sage: P5.Y(u1**-2*z)
            (-1.00000000000000)
            sage: z = UpperHalfPlaneProductElement([CC(1,3),CC(0,2)])
            sage: P5.Y(z)
            (0.210647934544181)
            sage: P5.Y(u1**2*z)
            (1.21064793454418)

            # The precision gets worse when the imaginary parts have large differences.
            sage: P5.Y(u1**2*z)
            (1.21064793454418)
            sage: P5.Y(u1**-2*z)
            (-0.789352065455819)
            sage: z = UpperHalfPlaneProductElement([CC(1,3),CC(0,2)])
            sage: P5.Y(z)
            (0.210647934544181)
            sage: P5.Y(u1**2*z)
            (1.21064793454418)
            sage: K4.<a> = NumberField(x^3-36*x-1)
            sage: H4 = ExtendedHilbertModularGroup(K4)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1, 1), CC(0, 2), CC(0, 1)])
            sage: P4.Y(z, return_error_estimate = True)
            ((-0.0644538365770865, 0.0000882959672881664), 2.220446049250313e-16)
            sage: H4 = ExtendedHilbertModularGroup(K4, tp_units = False)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.Y(z, return_error_estimate = True)
            ((-0.128907673154173, 0.0000882959672881664), 2.220446049250313e-16)
            sage: z = UpperHalfPlaneProductElement([CC(0, 100), CC(0, 200), CC(0, 100)])
            sage: H4 = ExtendedHilbertModularGroup(K4)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.Y(z, return_error_estimate = True)
            ((-0.0644538365770865, 0.0000882959672880350), 6.38378239159465e-16)
            sage: z = UpperHalfPlaneProductElement([CC(0, 1), CC(0, 2), CC(0, 100)])
            sage: P4.Y(z, return_error_estimate = True)
            ((0.319487796146047, 0.680877344221728), 1.199040866595169e-14)
            sage: CF = ComplexField(103)
            sage: z = UpperHalfPlaneProductElement([CF(0, 1), CF(0, 2), CF(0, 100)])
            sage: P4.Y(z, return_error_estimate = True)
            ((0.319487796146048307029200530746, 0.680877344221731788614078673856),
            1.2227344030925683e-29)
        """
        self._check_upper_half_plane_element(z)
        normalized_imag = z / z.imag_norm() ** (1 / self.group().number_field().degree())
        log_vector = matrix(vector(normalized_imag.imag_log())).transpose()
        B = self.basis_matrix_logarithmic_unit_lattice(prec=z.base_ring().prec())
        coordinate_vector = B.solve_right(log_vector, check=False)
        if return_error_estimate:
            return vector(coordinate_vector), (B * coordinate_vector - log_vector).norm(Infinity)
        return vector(coordinate_vector)

    def reduce_by_units(self, z, return_map=False):  # Doing
        r"""
        Reduce the point z with respect to action of units to a point where Im(z) belongs
        to a fundamental domain for the logarithmic unit lattice.

        INPUT:

        - ``z`` -- element of type UpperHalfPlaneProductElement
        - ``return_map`` -- boolean (default=False)
                            Set to True to return the map which does the reduction.

        EXAMPLES::
            sage: from hilbert_modgroup.all import *
            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: z = UpperHalfPlaneProductElement([CC(1, 1), CC(1, 1)]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: P5.reduce_by_units(z)
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: z = UpperHalfPlaneProductElement([CC(1, 1), CC(1, 1000)]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1000.00000000000*I]
            sage: P5.reduce_by_units(z, return_map = True)
            (
            [46.9787137637478 + 46.9787137637478*I, 0.0212862362522088 + 21.2862362522088*I],
            <BLANKLINE>
            [-21/2*a + 47/2              0]
            [             0              1]
            )
            sage: K3.<a> = QuadraticField(3)
            sage: lattice_ideal = K3.different()
            sage: level_ideal = K3.fractional_ideal(1)
            sage: H3 = ExtendedHilbertModularGroup(K3)
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: z = UpperHalfPlaneProductElement([CC(1, 0.01), CC(1, 0.002)]); z
            [1.00000000000000 + 0.0100000000000000*I, 1.00000000000000 + 0.00200000000000000*I]
            sage: P3.reduce_by_units(z, return_map = True)
            (
            [0.267949192431123 + 0.00267949192431123*I, 3.73205080756888 + 0.00746410161513775*I],
            <BLANKLINE>
            [a + 2     0]
            [    0     1]
            )

        """
        tp_units = self.group().tp_units()
        units = self.fundamental_units()
        # Only include the units != -1
        # To avoid overflow it is more efficient to apply the map,
        # e.g. compute (z*u**-k)/u**k instead of z*u**-(2k)
        if tp_units:
            floors = [-stable_floor(y) for y in self.Y(z)]
            reducing_map = prod([self.group().E(u**y) for u, y in zip(units, floors, strict=False)])
        else:
            floors = [-stable_floor(y / 2) for y in self.Y(z)]
            reducing_map = prod([self.group().E(u**y) for u, y in zip(units, floors, strict=False)])
        reduced_point = z.apply(reducing_map)
        if return_map:
            return reduced_point, reducing_map
        return reduced_point

    def is_reduced_by_units(self, z):
        r"""
        Return True if z is reduced with respect to the logarithmic unit lattice,
         otherwise return False.

        INPUT:

        - ``z`` -- element of the type UpperHalfPlaneProductElement_class

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,1)]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: P5.is_reduced_by_units(z)
            True
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(0,1000)]);z
            [1.00000000000000 + 1.00000000000000*I, 1000.00000000000*I]
            sage: P5.is_reduced_by_units(z)
            False
            sage: P5.reduce_by_units(z, return_map = True)
            (
            [46.9787137637478 + 46.9787137637478*I, 21.2862362522088*I],
            <BLANKLINE>
            [-21/2*a + 47/2              0]
            [             0              1]
            )
            sage: w,A=_
            sage: P5.is_reduced_by_units(w)
            True

        """
        tp_units = self.group().tp_units()
        if tp_units:
            return all(-1 / 2 <= y < 1 / 2 for y in self.Y(z))
        else:
            return all(-1 <= y < 1 for y in self.Y(z))

    def _construct_ideal(self, a, b=None):
        r"""
        Construct an ideal of the number field associated with self from a.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: P1._construct_ideal(1)
            Fractional ideal (1)
            sage: P1._construct_ideal(P1.number_field().fractional_ideal(1))
            Fractional ideal (1)
            sage: a = P1.number_field().fractional_ideal(1)
            sage: b = P1.number_field().fractional_ideal(2)
            sage: P1._construct_ideal(a, b)
            Fractional ideal (1/2)
        """
        if a is None:
            ideala = self.group().number_field().ring_of_integers().fractional_ideal(1)
        elif isinstance(a, NumberFieldIdeal):
            if b:
                ideala = a / b
            else:
                ideala = a
        elif a in self.group().number_field().ring_of_integers() and not b:
            ideala = self.group().number_field().ring_of_integers().fractional_ideal(a)
        elif (
            a in self.group().number_field().ring_of_integers()
            and b in self.group().number_field().ring_of_integers()
        ):
            ideala = self.group().ideal(a, b)
        else:
            raise ValueError(f"Could not construct a number field ideal from a={a} and b={b}")
        return ideala

    @cached_method
    def basis_matrix_ideal(self, a=None, prec=53):
        r"""
        Return the Basis matrix corresponding to an integer basis of an ideal a.

        INPUT:

        - ``a`` -- ideal or number field element.
        - ``prec`` -- integer (default=53)

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: P5.basis_matrix_ideal()
            [ 1.00000000000000 -1.61803398874989]
            [ 1.00000000000000 0.618033988749895]
            sage: P5.basis_matrix_ideal(2)
            [ 2.00000000000000 -3.23606797749979]
            [ 2.00000000000000  1.23606797749979]
            sage: K10.<a> = QuadraticField(10)
            sage: lattice_ideal = K10.different()
            sage: H10 = ExtendedHilbertModularGroup(K10)
            sage: P10 = ExtendedHilbertPullback(H10)
            sage: P10.basis_matrix_ideal()
            [ 1.00000000000000 -3.16227766016838]
            [ 1.00000000000000  3.16227766016838]

        """
        ideala = self._construct_ideal(a)
        entries = [list(beta.complex_embeddings(prec)) for beta in ideala.integral_basis()]
        n = self.group().number_field().degree()
        return matrix(RealField(prec), n, n, entries).transpose()

    # @cached_method
    # def basis_matrix_ideal_on_power_basis(self, a=None):
    #    r"""
    #    Return the Basis matrix corresponding to an integer basis of an ideal a
    #    in terms of the standard power basis.

    #    INPUT:

    #    - ``a`` -- ideal or number field element.

    #    EXAMPLES::

    #        sage: from hilbert_modgroup.extended.all import *
    #        sage: H1 = ExtendedHilbertModularGroup(5)
    #        sage: P1 = ExtendedHilbertPullback(H1)
    #        sage: P1.basis_matrix_ideal_on_power_basis()
    #        [   1 -1/2]
    #        [   0  1/2]
    #        sage: P1.basis_matrix_ideal_on_power_basis(2)
    #        [ 2 -1]
    #        [ 0  1]
    #        sage: H2=ExtendedHilbertModularGroup(10)
    #        sage: P2 = ExtendedHilbertPullback(H2)
    #        sage: P2.basis_matrix_ideal_on_power_basis()
    #        [1 0]
    #        [0 1]
    #        sage: a=H2.OK().gen(1)
    #        sage: P2.basis_matrix_ideal_on_power_basis(a+1)
    #        [9 1]
    #        [0 1]

    #    """
    #    ideala = self._construct_ideal(a)
    #    entries = [list(beta.vector()) for beta in ideala.integral_basis()]
    #    n = self.group().number_field().degree()
    #    return matrix(self.number_field(), n, n, entries).transpose()

    # def coordinates_in_number_field_ideal(self, x, a = None):
    #    r"""
    #    Return the coordinates of x with respect to an integral basis of a.

    #    INPUT:

    #    - ``x`` -- element of ideal a
    #    - ``a`` -- ideal or number field element.

    #    EXAMPLES::

    #        sage: from hilbert_modgroup.extended.all import *
    #        sage: H1 = ExtendedHilbertModularGroup(5)
    #        sage: P1 = ExtendedHilbertPullback(H1)
    #        sage: b1,b2=P1.number_field().fractional_ideal(1).basis()
    #        sage: P1.coordinates_in_number_field_ideal(b1)
    #        (1, 0)
    #        sage: P1.coordinates_in_number_field_ideal(b2)
    #        (0, 1)
    #        sage: H2 = ExtendedHilbertModularGroup(10)
    #        sage: P2 = ExtendedHilbertPullback(H2)
    #        sage: b1,b2 = P2.number_field().fractional_ideal(1).basis()
    #        sage: P2.coordinates_in_number_field_ideal(b1)
    #        (1, 0)
    #        sage: P2.coordinates_in_number_field_ideal(b2)
    #        (0, 1)

    #    """
    #    B = self.basis_matrix_ideal_on_power_basis(a = a)
    #    return B.inverse() * x.vector()

    # @cached_method
    # def basis_matrix_ideal__norm(self, a=None, prec=53, row=None):
    #    r"""
    #    Return the Basis matrix corresponding to an integer basis of an ideal a.

    #    INPUT:

    #    - ``a`` -- ideal or number field element.
    #    - ``prec`` -- integer (default=53)

    #    EXAMPLES::

    #        sage: from hilbert_modgroup.extended.all import *
    #        sage: H1 = ExtendedHilbertModularGroup(5)
    #        sage: P1 = ExtendedHilbertPullback(H1)
    #        sage: P1.basis_matrix_ideal__norm() # abs tol 1e-10
    #        2.61803398874989
    #        sage: P1.basis_matrix_ideal__norm(2) # abs tol 1e-10
    #        5.23606797749979
    #        sage: H2=ExtendedHilbertModularGroup(10)
    #        sage: P2 = ExtendedHilbertPullback(H2)
    #        sage: P2.basis_matrix_ideal__norm() # abs tol 1e-10
    #        4.16227766016838
    #        sage: a=H2.OK().gen(1)
    #        sage: P2.basis_matrix_ideal__norm(a+1) # abs tol 1e-10
    #        13.16227766016838

    #    """
    #    B = self.basis_matrix_ideal(a, prec=prec)
    #    if row is None:
    #        return B.norm(Infinity)
    #    elif 0 <= row < B.nrows():
    #        return sum([abs(x) for x in B.row(row)])
    #    else:
    #        raise ValueError(f"Can not find row:{row}")

    def X(self, z, a=None):
        r"""
        Coordinate of z with respect to the integral basis of an ideal.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.X(z)
            (1.00000000000000, 0.000000000000000)
            sage: b1,b2 = H1.base_ring().fractional_ideal(1).integral_basis(); b1,b2
            (1, 1/2*a - 1/2)
            sage: P1.X(z+b1)
            (2.00000000000000, 0.000000000000000)
            sage: P1.X(z-b1)
            (0.000000000000000, 0.000000000000000)
            sage: P1.X(z+b2)
            (0.999999999999999, 1.00000000000000)
            sage: P1.X(z-b2)
            (1.00000000000000, -1.00000000000000)

        """
        self._check_upper_half_plane_element(z)
        B = self.basis_matrix_ideal(a, prec=z.base_ring().prec())
        return vector(B**-1 * vector(z.real()))

    def reduce_by_translations(self, z, a=None, return_map=False):
        r"""
        Reduce the point z with respect to the cuspidal region in a neighbourhood of the cusp .

        INPUT:

        - ``z`` -- point in the upper half-plane.
        - ``a`` -- ideal cusp (default=None) if None use the ring of integers.
        - ``return_map`` -- boolean (default=False),
                            If set to True also return the map that makes the reduction.

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: from hilbert_modgroup.extended.all import *
            sage: K5.<a> = QuadraticField(5)
            sage: H5 = ExtendedHilbertModularGroup(K5)
            sage: P5 = ExtendedHilbertPullback(H5)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P5.reduce_by_translations(z, return_map=True)
            (
                                                      [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )
            sage: z = UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: P5.reduce_by_translations(z, return_map = True)
            (
            [-0.236067977499790 + 1.00000000000000*I, 0.236067977499790 + 1.00000000000000*I],
            <BLANKLINE>
            [    1 a - 1]
            [    0     1]
            )

        """
        X = self.X(z, a)
        for i in range(len(X)):
            if abs(abs(X[i]) - 0.5) < 1e-12:
                log.debug("CUSPIDAL WALL HIT: %d %s", i, X[i])
        ideala = self._construct_ideal(a)
        basis = ideala.integral_basis()
        correction = sum(b * stable_floor(X[i]) for i, b in enumerate(basis))
        reduced_point = z - correction
        if return_map:
            return reduced_point, Matrix(2, 2, [1, -correction, 0, 1])
        return reduced_point

    def is_reduced_by_translations(self, z, a=None, prec=53):
        r"""
        Check if the given point is reduced in the cuspidal region.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.is_reduced_by_translations(z)
            False
            sage: P1.reduce_by_translations(z)
            [1.00000000000000*I, 1.00000000000000*I]
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z))
            True
            sage: z = UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: P1.is_reduced_by_translations(z)
            False
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z))
            True
            sage: b1,b2 = H1.OK().gens(); b1+b2
            3/2*a + 1/2
            sage: P1.is_reduced_by_translations(z,b1+b2)
            False
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z,b1+b2),b1+b2)
            True
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4.is_reduced_by_translations(z) # abs tol 1e-10
            False
            sage: P4.is_reduced_by_translations(P4.reduce_by_translations(z))
            False

        """
        X = self.X(z, a)
        return all(-1 / 2 <= x < 1 / 2 for x in X)

    def basis_matrix_ideal_plusz(self, z, a=None):
        r"""
        Return the Basis matrix corresponding to the lattice
        L = OKz + OK embedded in R^{2n} and given by the fixed
        integral basis of the ideal a.

        INPUT:

        - ``z`` -- element of a product of complex planes
        - ``a`` -- ideal or number field element.
        - ``prec`` -- integer (default=53)

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1.basis_matrix_ideal_plusz(z, P1._construct_ideal(1))
            [  1.00000000000000   1.00000000000000  0.000000000000000  0.000000000000000]
            [ -1.61803398874989  0.618033988749895  0.000000000000000  0.000000000000000]
            [ 0.000000000000000  0.000000000000000   1.00000000000000   1.00000000000000]
            [-0.000000000000000  0.000000000000000  -1.61803398874989  0.618033988749895]
            sage: z = UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
            sage: P1.basis_matrix_ideal_plusz(z,P1._construct_ideal(1))
            [  1.00000000000000   1.00000000000000  0.000000000000000  0.000000000000000]
            [ -1.61803398874989  0.618033988749895  0.000000000000000  0.000000000000000]
            [ 0.000000000000000  0.000000000000000  0.500000000000000   1.00000000000000]
            [-0.000000000000000  0.000000000000000 -0.809016994374947  0.618033988749895]
            sage: z = UpperHalfPlaneProductElement([CC(2.2,0.5),CC(1,0.5)])
            sage: P1.basis_matrix_ideal_plusz(z,P1._construct_ideal(1))
            [  1.00000000000000   1.00000000000000  0.000000000000000  0.000000000000000]
            [ -1.61803398874989  0.618033988749895  0.000000000000000  0.000000000000000]
            [  2.20000000000000   1.00000000000000  0.500000000000000  0.500000000000000]
            [ -3.55967477524977  0.618033988749895 -0.809016994374947  0.309016994374947]
            sage: K2.<a> = QuadraticField(10)
            sage: H2 = ExtendedHilbertModularGroup(K2)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(2,0.5),CC(1,0.2)])
            sage: P2.basis_matrix_ideal_plusz(z,P2._construct_ideal(1))
            [ 1.00000000000000  1.00000000000000 0.000000000000000 0.000000000000000]
            [-3.16227766016838  3.16227766016838 0.000000000000000 0.000000000000000]
            [ 2.00000000000000  1.00000000000000 0.500000000000000 0.200000000000000]
            [-6.32455532033676  3.16227766016838 -1.58113883008419 0.632455532033676]
            sage: zv = [23.3400000000000 + 0.0100000000000000*I,
            ....:       0.0200000000000000 + 0.0300000000000000*I]
            sage: z = UpperHalfPlaneProductElement(zv)
            sage: P2.basis_matrix_ideal_plusz(z)
            [   1.00000000000000    1.00000000000000   0.000000000000000   0.000000000000000]
            [  -3.16227766016838    3.16227766016838   0.000000000000000   0.000000000000000]
            [   23.3400000000000  0.0200000000000000  0.0100000000000000  0.0300000000000000]
            [  -73.8075605883300  0.0632455532033676 -0.0316227766016838  0.0948683298050514]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4.basis_matrix_ideal_plusz(z)
            [ 1.00000000000000  1.00000000000000 0.000000000000000 0.000000000000000]
            [-1.73205080756888  1.73205080756888 0.000000000000000 0.000000000000000]
            [ 1.00000000000000  2.00000000000000  1.00000000000000  1.00000000000000]
            [-1.73205080756888  3.46410161513775 -1.73205080756888  1.73205080756888]


        """
        ideala = self._construct_ideal(a)
        zero = self.number_field()(0)
        prec = z.prec()
        # In case it is an Upper-half plane element we need to allow multiplication
        # by elements on the real axis.
        z = z.as_ComplexPlaneProductElement()
        entries = [
            beta.complex_embeddings(prec) + zero.complex_embeddings(prec)
            for beta in ideala.integral_basis()
        ]
        entries += [
            (z * z.parent()(beta)).real() + (z * z.parent()(beta)).imag()
            for beta in ideala.integral_basis()
        ]
        n = self.group().base_ring().degree()
        return matrix(RealField(prec), 2 * n, 2 * n, entries)

    def _shortest_vectors_ideal_plusz(self, z, a=None, return_scaled_matrix=False):
        r"""
        Compute a list of potentially shortest vectors in the lattice az+a using LLL.

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``a`` -- ideal (default =None) if None we use the entire ring of integers as lattice.
        - ``return_scaled_matrix`` -- boolean (default False)
                                      Set to True to return the scaled matrix.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: K.<a> = QuadraticField(5)
            sage: H1 = ExtendedHilbertModularGroup(K)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1._shortest_vectors_ideal_plusz(z)
            [1 0 0 0]
            [0 0 1 0]
            [0 1 0 0]
            [0 0 0 1]
            sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
            sage: P1._shortest_vectors_ideal_plusz(z)
            [0 0 1 0]
            [0 0 0 1]
            [1 0 0 0]
            [1 1 0 0]
            sage: z=UpperHalfPlaneProductElement([CC(2.2,0.5),CC(1,0.5)])
            sage: P1._shortest_vectors_ideal_plusz(z,P1._construct_ideal(1))
            [-1  1  1  0]
            [ 0 -2  0  1]
            [ 1 -2  0  1]
            [ 2  0 -1  0]
            sage: K.<a> = QuadraticField(10)
            sage: H2 = ExtendedHilbertModularGroup(K)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(2,0.5),CC(1,0.2)])
            sage: P2._shortest_vectors_ideal_plusz(z,P2._construct_ideal(1))
            [-1  0  1  0]
            [ 2  0 -1  0]
            [-3 -1  3  1]
            [-5  2  2 -1]
            sage: z=UpperHalfPlaneProductElement([23.3400000000000 + 0.0100000000000000*I,
            ....:                                  0.0200000000000000 + 0.0300000000000000*I])
            sage: P2._shortest_vectors_ideal_plusz(z)
            [  25   -8    1    1]
            [  10   -3   -4   -1]
            [ 212  -67  -15    1]
            [ 348 -110  -14    5]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H3 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P3._shortest_vectors_ideal_plusz(z)
            [ 1  0  0  0]
            [-2  0  1  0]
            [ 2  1 -1  0]
            [-1 -1  0  1]


        NOTE:
            The idea behind this function was taken from [BoSt2015] (especially the scaling factor)

        REFERENCES:
            [BoSt2015] F. Boyer and M. Streng,
                       "Examples of CM curves of genus two defined over the reflex field",
                       LMS Journal of Comp. Math., Vol. 18 (2015), issue 01, pp 507-538\n",



        """
        basis_matrix = self.basis_matrix_ideal_plusz(z, a)
        n = basis_matrix.nrows() // 2
        # This is essentially the inverse of "machine epsilon"
        epsilon_inverse = 2 ** (z[0].parent().prec())
        # Scale all rows by epsilon_inverse*y^{-1/2} (we assume that N(y) is small)
        scaled_rows = []
        for basis_row in basis_matrix:
            row = []
            for i, b in enumerate(basis_row):
                entry = (epsilon_inverse * b / (z[i % n].imag().sqrt())).round()
                # entry = (epsilon_inverse * b).round()
                row.append(entry)
            scaled_rows.append(row)
        integral_scaled_basis_matrix = Matrix(ZZ, scaled_rows)
        if return_scaled_matrix:
            return integral_scaled_basis_matrix
        # Apply LLL to find a reduced basis
        R, U = integral_scaled_basis_matrix.LLL(transformation=True)
        # The shortest basis vector should be the first but in practice
        # one of the other vectors can provide a better estimate.
        return U

    def get_heuristic_closest_cusp(self, z, a=None, as_cusp=False):
        """
        Try to find a heuristic closest cusp using LLL.

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``a`` -- ideal or number field element (default = None)
          If None then this is set to the entire ring of integers.
        - ``as_cusp`` -- boolean (default: False)
          If True we return an element of type NFCusp_wrt_lattice_ideal

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: K.<a> = QuadraticField(5)
            sage: H1 = ExtendedHilbertModularGroup(K)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1.get_heuristic_closest_cusp(z)
            (-1, 0)
            sage: z = UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
            sage: P1.get_heuristic_closest_cusp(z)
            (0, 1)
            sage: z = UpperHalfPlaneProductElement([CC(2.2,0.5),CC(1,0.5)])
            sage: P1.get_heuristic_closest_cusp(z)
            (-1/2*a + 3/2, 1)
            sage: K.<a> = QuadraticField(10)
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z = UpperHalfPlaneProductElement([CC(2,0.5),CC(1,0.2)])
            sage: P2.get_heuristic_closest_cusp(z)
            (1, 1)
            sage: z = UpperHalfPlaneProductElement([23.3400000000000 + 0.0100000000000000*I,
            ....:                                   0.0200000000000000 + 0.0300000000000000*I])
            sage: P2.get_heuristic_closest_cusp(z)
            (3*a - 10, -a - 4)
            sage: w=UpperHalfPlaneProductElement([-0.668903800800698 + 0.0362571615120737*I,
            ....:                                  0.708560139622790 + 0.00414937759336099*I])
            sage: P2.get_heuristic_closest_cusp(w)
            (-a + 4, 2*a - 5)
            sage: x = ZZ['x'].gen()
            sage: H3=ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=ExtendedHilbertPullback(H3)
            sage: K3=P3.number_field()
            sage: z = UpperHalfPlaneProductElement([0+0.02*I,10+0.2*I,1+0.2*I])
            sage: P3.get_heuristic_closest_cusp(z)
            (1/3*a^2 + 1/3*a - 29/3, a - 1)
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4.get_heuristic_closest_cusp(z)
            (-1, 0)


        """
        a = self._construct_ideal(a)
        try:
            shortest_basis_vectors = self._shortest_vectors_ideal_plusz(z, a)
        except ValueError:
            log.critical(
                "The LLL 'finding shortest vector' has failed. "
                + "It is likely that you need to upgrade your version of Sage to 9.4+."
            )
            return None
        # Convert this lattice vector to two integers sigma and rho:
        n = len(a.basis())
        dmin = None
        cusp_min = None
        for vec in shortest_basis_vectors:
            rho = sum([vec[i] * a.basis()[i] for i in range(n)])
            sigma = sum([vec[n + i] * a.basis()[i] for i in range(n)])
            # The actual cusp that minimizes |sigma*z+rho| is of course -rho/sigma
            d = distance_to_cusp_eg(self, -rho, sigma, z)
            if not dmin or d < dmin:
                dmin = d
                cusp_min = -rho, sigma
        if as_cusp:
            cusp_min = self._construct_cusp(cusp_min)
        return cusp_min

    def _construct_cusp(self, c, d=None):
        r"""
        Return an instance of NFCusp_wrt_base_field for the number field of self from input c and optional d

        INPUT:

        - ``c`` -- instance of NFCusp_wrt_lattice_ideal, element or tuple of elements of number field
        - ``d`` -- element of number field or None (default: None)

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: K.<a> = QuadraticField(3)
            sage: P1 = HilbertPullback(HilbertModularGroup(K))
            sage: P1._construct_cusp(1,0)
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
            sage: P1._construct_cusp((1,0))
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
            sage: P1._construct_cusp(0, 1)
            Cusp [0: 1] of Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
        """
        lattice_ideal = self.group().lattice_ideal()
        if isinstance(c, NFCusp_wrt_lattice_ideal) and c.number_field() == self.number_field():
            return c
        if isinstance(c, NFCusp_wrt_lattice_ideal) and c.number_field() != self.number_field():
            raise ValueError(f"The input cusp {c} has wrong base number field.")
        if isinstance(c, tuple) and len(c) == 2:
            c, d = c
        try:
            cusp = NFCusp_wrt_lattice_ideal(lattice_ideal, c, d)
        except Exception as err:
            raise ValueError(
                f"Could not construct a cusp with respect to the lattice_ideal = {lattice_ideal} from c={c} and d={d}"
            ) from err
        return cusp

    def number_field(self):
        """
        Return number field of the Extended Hilbert modular group of self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: H3 = ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: P3.number_field()
            Number Field in a with defining polynomial x^3 - 36*x - 1

        """
        return self.group().number_field()

    def reduce_in_cuspidal_region(self, z, cusp=None, return_map=False):  # needcheking
        r"""
        Reduce the point z with respect to the cuspidal region in a neighbourhood of a
        representative cusp.

        INPUT:

        - ``z`` -- point in the upper half-plane.
        - ``cusp`` -- cusp of the group of self.
        - ``return_map`` -- boolean (default=False)
          Return the map A such that AN^-1z is reduced where N is the
          cusp normalizing map for the cusp.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.reduce_in_cuspidal_region(z, H1.cusps()[0])
            [1.00000000000000*I, 1.00000000000000*I]
            sage: P1.reduce_in_cuspidal_region(z,H1.cusps()[0],return_map=True)
            (
                                                      [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )
            sage: P1.reduce_in_cuspidal_region(z,P1._construct_cusp(0,1),return_map=True)
             (
                                                      [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )

        Check that if we apply a cusp-normalizing map to a point then the reduction with respect
        to that cusp is the same as the reduction of the original point with respect to infinity:
            sage: c = NFCusp_wrt_lattice_ideal(P1.group().lattice_ideal(), 1, 2)
            sage: N = H1.cusp_normalizing_map(c)
            sage: w = z.apply(N); w
            [0.384615384615385 + 0.0769230769230769*I, 0.384615384615385 + 0.0769230769230769*I]
            sage: P1.reduce_in_cuspidal_region(w, cusp = c,return_map = True)
            (
            [0.384615384615385 + 0.0769230769230769*I, 0.384615384615385 + 0.0769230769230769*I],
            [1 0]
            [0 1]
            )
            sage: P1.reduce_in_cuspidal_region(z, cusp = c, return_map = True)
             (
                                                       [ 1 -1]
             [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
             )

            sage: z = UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: P1.reduce_in_cuspidal_region(z)
            [-0.236067977499790 + 1.00000000000000*I, 0.236067977499790 + 1.00000000000000*I]
            sage: P1.reduce_in_cuspidal_region(z, P1.group().cusps()[0], return_map=True)
            (
            [-0.236067977499790 + 1.00000000000000*I, 0.236067977499790 + 1.00000000000000*I],
            <BLANKLINE>
            [     1 a - 1]
            [     0      1]
            )
            sage: H3 = ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1),CC(10,1)])
            sage: cusp = P3.group().cusps()[0]
            sage: P3.reduce_in_cuspidal_region(z, cusp,  # abs tol 1e-10
            ....:                   return_map=True)
            (
            [2.98606258583498 + 1.0*I, -1.97222162680976 + 1.0*I, -0.0138409590252291 + 1.0*I],
            <BLANKLINE>
            [     1 -a - 4]
            [     0      1]
            )
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: cusp = P4.group().cusps()[0]
            sage: P4.reduce_in_cuspidal_region(z, cusp, return_map = True)
            (
                                                                         [ 1 -1]
            [1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I], [ 0  1]
            )

        """
        self._check_upper_half_plane_element(z)
        if cusp is None:
            cusp = self.group().cusps()[0]
        ideala = ((self.group().lattice_ideal()) ** -1) * (cusp.ideal() ** -2)
        # Then reduce with respect to the units, followed by reduction by translation with respect
        # to the ideal (lattice_ideal**-1)*(a**-2)
        if return_map:
            w, A = self.reduce_by_units(z, return_map=return_map)
            w, B = self.reduce_by_translations(w, ideala, return_map=return_map)
            return w, B * A
        else:
            w = self.reduce_by_units(z, return_map=return_map)
            w = self.reduce_by_translations(w, ideala, return_map=return_map)
        return w

    def find_closest_cusp(self, z, return_multiple=False, as_cusp=True):  # Doing
        r"""
        Find the closest cusp (rho:sigma) to z

        INPUT:

        - `z` -- point in the upper half-plane
        - `return_multiple` -- boolean: default False
                               Set to True to return all cusps with the same minimal distance,
                                        otherwise just return one closest cusp.
        - ``as_cusp`` -- boolean (default True) return instance(s) of NFcusps or tuple(s)

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.find_closest_cusp(z, as_cusp=False)
            (1, 0)
            sage: P1.find_closest_cusp(z)
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? with respect to lattice_ideal
            sage: P1.find_closest_cusp(z, return_multiple=True)
            [Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
             with respect to lattice_ideal, Cusp [-1: -1] of Number Field in a with defining polynomial
             x^2 - 5 with a = 2.236067977499790? with respect to  lattice_ideal]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: K2 = H2.OK().number_field()
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(0,1.0),CC(0,1.0)])
            sage: P2.find_closest_cusp(z)
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380?
            with respect to lattice_ideal
            sage: P2.find_closest_cusp(z,return_multiple=True)
            [Cusp Infinity of Number Field in a with defining polynomial x^2 - 10 with
            a = 3.162277660168380? with respect to lattice_ideal, Cusp [0: 1] of Number Field in a with
            defining polynomial x^2 - 10 with a = 3.162277660168380? with respect to  lattice_ideal]
            sage: z=UpperHalfPlaneProductElement([CC(2.58,0.5),CC(0.5,0.5)])
            sage: P2.find_closest_cusp(z)
            Cusp [-a: -a - 2] of Number Field in a with defining polynomial x^2 - 10 with
            a = 3.162277660168380? with respect to  lattice_ideal
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4.find_closest_cusp(z)
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 3 with a = 1.732050807568878?
            with respect to lattice_ideal

        """
        closest_cusp = find_closest_cusp(
            self, z, return_multiple=return_multiple, use_lll=True, use_norm_bound=True
        )
        lattice_ideal = self.group().lattice_ideal()
        if as_cusp and return_multiple:
            return [NFCusp_wrt_lattice_ideal(lattice_ideal, c[0], c[1]) for c in closest_cusp]
        if as_cusp:
            return NFCusp_wrt_lattice_ideal(lattice_ideal, closest_cusp[0], closest_cusp[1])
        return closest_cusp

    def distance_to_cusp_eg(self, cusp, z):
        """
        Give the distance from the point z to the cusp: N(a)^-1 N(Im A^-1z)^(-1/2)

        INPUT:
        - `cusp` -- NF cusp wrt lattice_ideal or tuple of integral elements in the number field.
        - `z`

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = ComplexPlaneProductElement([CC(0,1),CC(0,1)]);
            sage: P1.distance_to_cusp_eg(H1.cusps()[0], z)
            1.00000000000000
            sage: P1.distance_to_cusp_eg(H1.cusps()[0], z*[2,2])
            0.500000000000000
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z =  ComplexPlaneProductElement([CC(0,1),CC(0,1)]);
            sage: P2.distance_to_cusp_eg(H2.cusps()[0], z)
            1.00000000000000
            sage: P2.distance_to_cusp_eg(H2.cusps()[1], z)
            62.7694193059008
            sage: c = NFCusp_wrt_lattice_ideal(H2.lattice_ideal(), 3, H2.lattice_ideal().number_field().gen()+1)
            sage: P2.distance_to_cusp_eg(c, z) # abs tol 1e-10
            6.32455532033676
            sage: z1 = ComplexPlaneProductElement([CC(-1.38,0.1),CC(0.72,0.1)]); z1
            [-1.38000000000000 + 0.100000000000000*I, 0.720000000000000 + 0.100000000000000*I]
            sage: P2.distance_to_cusp_eg(H2.cusps()[0], z1)
            10.0000000000000
            sage: c=NFCusp_wrt_lattice_ideal(H2.lattice_ideal(), 3,H2.lattice_ideal().number_field().gen()+1)
            sage: P2.distance_to_cusp_eg(c, z1) # abs tol 1e-10
            0.300834689631305
            sage: P2.distance_to_cusp_eg(H2.cusps()[1], z1)
            831.857032071542
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: c = NFCusp_wrt_lattice_ideal(H4.lattice_ideal(), 3, H4.lattice_ideal().number_field().gen()+1)
            sage: P4.distance_to_cusp_eg(H4.cusps()[0], z1)
            10.0000000000000

        """
        cusp = self._construct_cusp(cusp)
        if isinstance(z, ComplexPlaneProductElement__class):
            z = UpperHalfPlaneProductElement(z.z())
        elif not isinstance(z, UpperHalfPlaneProductElement__class):
            z = UpperHalfPlaneProductElement(z)  # Try to make an upper half-plane element
        return distance_to_cusp_eg(self, cusp.numerator(), cusp.denominator(), z)

    def polytope_from_bounds(self, bounds, B=None):
        r"""
        Return the polytope defined by a set of bounds, either after applying the map 'B', or not.

        Get a list of integer points in an ideal with embeddings in RR^n within a given bound.
        Reference: Algorithm 10

        INPUT:

        - ``a`` -- ideal or algebraic integer.
        - ``bounds`` -- list of bounds for the coordinates or scalar which then gives a cube with
                        the same bound.
                        Can be of the form [b1,b2,...,bn] or [(a1,b1),(a2,b2)....(an,bn)].
                        In the first case the bounds are interpreted as (-b1,b1),...,(-bn,bn).
        - ``return_polyhedron`` -- boolean (default False)
                                   Set to True to return a polyhedron of the corresponding domain
        - ``preimage`` -- boolean (default False)
                          Set to True to return the polyhedron of the pre-image
                          (only used when return_polyhedron=True)

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: p = P1.polytope_from_bounds([(-1,1),(-1,1)])
            sage: p.vertices()
            (A vertex at (-1.0, -1.0),
             A vertex at (-1.0, 1.0),
             A vertex at (1.0, -1.0),
             A vertex at (1.0, 1.0))
            sage: p = P1.polytope_from_bounds([1,1])
            sage: p.vertices()
            (A vertex at (-1.0, -1.0),
             A vertex at (-1.0, 1.0),
             A vertex at (1.0, -1.0),
             A vertex at (1.0, 1.0))

        """
        from sage.all import Polyhedron, vector

        if not isinstance(bounds, list | tuple):
            raise ValueError("Need a list of bounds!")
        if not isinstance(bounds[0], tuple):
            bounds = [(-b, b) for b in bounds]
        # Hypercube we want embeddings in
        vertices = cartesian_product([[RDF(a), RDF(b)] for a, b in bounds]).list()
        p1 = Polyhedron(vertices, base_ring=RDF)
        if not B:
            return p1
        # Hypercube containing the integral points of the coordinates wrt the integral basis
        vertices = [vector(list(B * vector(x))) for x in p1.vertices()]
        # Try to make a polyhedron of the mapped vertices.
        return Polyhedron(vertices, base_ring=RDF)

    @cached_method
    def max_ideal_norm(self):  # Need to work
        r"""
        Compute the maximum of the norms of all ideal representatives in the ideal class group.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: P1.max_ideal_norm()
            1
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: P2.max_ideal_norm()
            2
            sage: H3=ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=ExtendedHilbertPullback(H3)
            sage: P3.max_ideal_norm()
            6
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.max_ideal_norm()
            1

        """
        return max([x.norm() for x in self.group().ideal_cusp_representatives()])

    def _matrix_BLambda_row_sum(self, i=None):  # Need to work
        r"""
        Compute ``r_i(B_{\Lambda}) = sum_j |b_ij|`` or ``sum_ij |b_ij|``.

        INPUT:

        -`` i`` -- integer (default: None)
                    If i is given then return sum of absolute values in row nr. i,
                    otherwise return the sum.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: P1._matrix_BLambda_row_sum() # abs tol 1e-10
            0.962423650119207
            sage: P1._matrix_BLambda_row_sum(0) # abs tol 1e-10
            0.481211825059603
            sage: P1._matrix_BLambda_row_sum(1) # abs tol 1e-10
            0.481211825059603

            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: P2._matrix_BLambda_row_sum() # abs tol 1e-10
            3.63689291846413
            sage: P2._matrix_BLambda_row_sum(0) # abs tol 1e-10
            1.81844645923207
            sage: P2._matrix_BLambda_row_sum(1) # abs tol 1e-10
            1.81844645923207

            sage: H3 = ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: P3._matrix_BLambda_row_sum()
            15.7133517843634
            sage: P3._matrix_BLambda_row_sum(0)
            6.06261225308888
            sage: P3._matrix_BLambda_row_sum(1)
            5.37061649381160
            sage: P3._matrix_BLambda_row_sum(2)
            4.28012303746290
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4._matrix_BLambda_row_sum()
            2.63391579384963
            sage: P4._matrix_BLambda_row_sum(0)
            1.31695789692482
            sage: P4._matrix_BLambda_row_sum(1)
            1.31695789692482

        """
        K = self.number_field()
        lattice_ideal = self.group().lattice_ideal()
        # level_ideal = self.group().level_ideal()
        H = ExtendedHilbertModularGroup(K, lattice_ideal=lattice_ideal, tp_units=False)
        P = ExtendedHilbertPullback(H)
        B = P.basis_matrix_logarithmic_unit_lattice()
        # tp_units = self.group().tp_units()
        # if tp_units:
        #    if i is not None:
        #        return 1/2* sum([abs(x) for x in B[i]])
        #    else:
        #        return 1/2 * sum([sum([abs(x) for x in row]) for row in B])
        # else:
        if i is not None:
            return sum([abs(x) for x in B[i]])
        else:
            return sum([sum([abs(x) for x in row]) for row in B])

    @cached_method
    def _exp_matrix_BLambda_row_sum(self, i=None):  # Doing
        r"""
        Compute ``exp(r_i(B_{\Lambda})) = sum_j |b_ij|`` or ``sum_ij |b_ij|``.

        INPUT:

        -`` i`` -- integer (default: None)
                   If i is given then return sum of absolute values in row nr. i,
                    otherwise return the sum.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: P1._exp_matrix_BLambda_row_sum()
            2.61803398874989
            sage: P1._exp_matrix_BLambda_row_sum(0)
            1.61803398874989
            sage: P1._exp_matrix_BLambda_row_sum(1)
            1.61803398874989

            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: P2._exp_matrix_BLambda_row_sum() # abs tol 1e-10
            37.9736659610102
            sage: P2._exp_matrix_BLambda_row_sum(0) # abs tol 1e-10
            6.16227766016838
            sage: P2._exp_matrix_BLambda_row_sum(1) # abs tol 1e-10
            6.16227766016837

            sage: H3=ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=ExtendedHilbertPullback(H3)
            sage: P3._exp_matrix_BLambda_row_sum() # abs tol 1e-7
                6.67147667780288e6
            sage: P3._exp_matrix_BLambda_row_sum(0) # abs tol 1e-10
                429.495924779268
            sage: P3._exp_matrix_BLambda_row_sum(1) # abs tol 1e-10
                214.995370171897
            sage: P3._exp_matrix_BLambda_row_sum(2) # abs tol 1e-10
                72.2493288346009
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4._exp_matrix_BLambda_row_sum()
            13.9282032302755
            sage: P4._exp_matrix_BLambda_row_sum(0)
            3.73205080756888
            sage: P4._exp_matrix_BLambda_row_sum(1)
            3.73205080756888


        """
        return self._matrix_BLambda_row_sum(i).exp()

    @cached_method()
    def Di(self, i=None):  # Doing
        r"""
        Return the bound ``D_i`` for this Hilbert modular group.

        INPUT:

        -`` i`` -- integer (default: None)
                    If i is given then return the bound with exp(r_i(B_Lambda))
                    else use exp(sum_i r_i(B_Lambda)

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: from hilbert_modgroup.extended.all import *
            sage: P1 = HilbertPullback(HilbertModularGroup(5))
            sage: P1.Di() # abs tol 1e-10
            1.6180339887498947
            sage: P1.Di(0) # abs tol 1e-10
            1.272019649514069
            sage: P1.Di(1) # abs tol 1e-10
            1.272019649514069
            sage: P2 = HilbertPullback(HilbertModularGroup(10))
            sage: P2.Di() # abs tol 1e-10
            8.714776642118862
            sage: P2.Di(0) # abs tol 1e-10
            3.510634603648856
            sage: P2.Di(1) # abs tol 1e-10
            3.510634603648854
            sage: x = ZZ['x'].gen()
            sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
            sage: P3 = HilbertPullback(HilbertModularGroup(K3))
            sage: P3.Di() # abs tol 1e-10
            5.048917339522305
            sage: P3.Di(0) # abs tol 1e-10
                1.4989928631309308
            sage: P3.Di(1) # abs tol 1e-10
                1.6738989622449851
            sage: P3.Di(2) # abs tol 1e-10
                2.0121921726123237
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.Di()
            3.7320508075688767
            sage: P4.Di(0)
            1.9318516525781362
            sage: P4.Di(1)
            1.9318516525781364
        """
        n = self.group().number_field().degree()
        return float(self.max_ideal_norm() ** (1 / n) * self._exp_matrix_BLambda_row_sum(i).sqrt())

    @cached_method()
    def D(self):
        r"""
        Return a list of all bounds ``D_i`` for this Hilbert modular group.


        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: P1 = ExtendedHilbertPullback(ExtendedHilbertModularGroup(5))
            sage: P1.D() # abs tol 1e-10
            [1.272019649514069, 1.272019649514069]
            sage: P2 = ExtendedHilbertPullback(ExtendedHilbertModularGroup(10))
            sage: P2.D() # abs tol 1e-10
            [3.510634603648856, 3.510634603648854]
            sage: x = ZZ['x'].gen()
            sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
            sage: P3 = ExtendedHilbertPullback(ExtendedHilbertModularGroup(K3))
            sage: P3.D() # abs tol 1e-10
                [1.4989928631309308, 1.6738989622449851, 2.0121921726123237]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4.D()
            [1.9318516525781362, 1.9318516525781364]
        """
        n = self.group().number_field().degree()
        return [self.Di(i) for i in range(n)]

    def _bound_for_closest_cusp(self):
        """
        This is the bound such that if a cusp is closer than this then it is the closest cusp.
        Reference:

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: P1 = ExtendedHilbertPullback(ExtendedHilbertModularGroup(5))
            sage: P1._bound_for_closest_cusp() # abs tol 1e-10
            0.190983005625053
            sage: P2 = ExtendedHilbertPullback(ExtendedHilbertModularGroup(10))
            sage: P2._bound_for_closest_cusp() # abs tol 1e-10
            0.00658350974743101
            sage: x = ZZ['x'].gen()
            sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
            sage: P3 = ExtendedHilbertPullback(ExtendedHilbertModularGroup(K3))
            sage: P3._bound_for_closest_cusp() # abs tol 1e-10
            0.0138694259275406
            sage: K.<a> = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4._bound_for_closest_cusp()
            0.0358983848622454

        """
        n = self.group().number_field().degree()
        return self.max_ideal_norm() ** (-1) * 2 ** (-n / 2.0) / self._exp_matrix_BLambda_row_sum()

    def _Dzi(self, z, i, initial_bd_d=None, use_initial_bd_d=True):
        """
        Return the constant `a_i` used to bound the embeddings of sigma.

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``initial_bd_d`` -- an initial bound for the distance to nearest cusp (default None)
        - ``use_initial_bd_d`` -- boolean (default: 'True') Use the initial bound or not.
                                This should only be set to False for demonstration or testing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1._Dzi(z,0)
            1.27201964951407
            sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
            sage: P1._Dzi(z,0)
            1.06963676340867
            sage: z=UpperHalfPlaneProductElement([CC(2.2,0.5),CC(1,0.5)])
            sage: P1._Dzi(z,0)
            1.79890743994787
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(2,0.5),CC(1,0.2)])
            sage: P2._Dzi(z,0)
            6.24288923183891
            sage: z=UpperHalfPlaneProductElement([23.3400000000000 + 0.0100000000000000*I,
            ....:                                  0.0200000000000000 + 0.0300000000000000*I])
            sage: P2._Dzi(z,0)
            24.4704299162552
            sage: x = ZZ['x'].gen()
            sage: H3=ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=ExtendedHilbertPullback(H3)
            sage: K3=P3.number_field()
            sage: z = UpperHalfPlaneProductElement([0+0.02*I,10+0.2*I,1+0.2*I])
            sage: P3._Dzi(z,0)
            72.7600147269406
            sage: P3._Dzi(z,1)
            51.4787281347537
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4._Dzi(z, 0)
            1.93185165257814

        """
        n = self.group().number_field().degree()
        if not use_initial_bd_d:
            return self.Di(i) * z.imag_norm() ** (-1 / (2 * n))
        dist_to_infinity_bd = z.imag_norm() ** (-1 / 2)
        dist_to_zero_bd = ((self.group().lattice_ideal().inverse().norm()) ** -1) * (
            z.abs_square_norm() / z.imag_norm()
        ) ** (0.5)  # Need to work
        if initial_bd_d:
            d = min(initial_bd_d, dist_to_infinity_bd, dist_to_zero_bd)
        else:
            d = min(dist_to_infinity_bd, dist_to_zero_bd)
        d = d ** (1 / n)
        return self.Di(i) * d

    def _Dz(self, z, initial_bd_d=None, use_initial_bd_d=True):
        """
        Return the vector of all bounds Dzi

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``initial_bd_d`` -- an initial bound for the distance to nearest cusp (default None)
        - ``use_initial_bd_d`` -- boolean (default: 'True') Use the initial bound or not.
                                This should only be set to False for demonstration or testing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1._Dz(z)
            [1.27201964951407, 1.27201964951407]
            sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
            sage: P1._Dz(z)
            [1.06963676340867, 1.06963676340867]
            sage: z=UpperHalfPlaneProductElement([CC(2.2,0.5),CC(1,0.5)])
            sage: P1._Dz(z)
            [1.79890743994787, 1.79890743994787]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(2,0.5),CC(1,0.2)])
            sage: P2._Dz(z) # abs tol 1e-10
            [6.24288923183892, 6.24288923183892]
            sage: z=UpperHalfPlaneProductElement([23.3400000000000 + 0.0100000000000000*I,
            ....:                                  0.0200000000000000 + 0.0300000000000000*I])
            sage: P2._Dz(z) # abs tol 1e-10
            [24.4704299162553, 24.4704299162552]
            sage: x = ZZ['x'].gen()
            sage: H3 = ExtendedHilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3 = ExtendedHilbertPullback(H3)
            sage: K3= P3.number_field()
            sage: z = UpperHalfPlaneProductElement([0+0.02*I,10+0.2*I,1+0.2*I])
            sage: P3._Dz(z)
            [72.7600147269406, 51.4787281347537, 29.8421537172252]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4._Dz(z)
            [1.93185165257814, 1.93185165257814]
        """
        n = self.number_field().degree()
        return [
            self._Dzi(z, i, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d)
            for i in range(n)
        ]

    @cached_method()
    def _bound_for_sigma_norm(self, z, dist=None):
        r"""
        Return the bound for the norm of sigma.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_sigma_norm(z) # abs tol 1e-8
            1.1111111111111112
            sage: P1._bound_for_sigma_norm(z, 1) # abs tol 1e-8
            1.0540925533894598
        """
        if not dist:
            d = z.imag_norm() ** (-1)
        else:
            d = dist * z.imag_norm() ** (-1 / 2)
        return float(d * self.max_ideal_norm())

    def _bound_for_sigma_embeddings(self, z, initial_bd_d=None, prec=16, use_initial_bd_d=True):
        """
        Bound for the embeddings of the denominator of the closest cusp to z


        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``prec`` -- integer - the number of bits precision in the returned values.
        - ``initial_bd_d`` -- float - an initial bound (default None) for the distance to a cusp.
                              If it is None or larger than the distance to infinity then this
                              distance is used.
        - ``use_initial_bd_d`` -- boolean (default: 'True') -- Use the initial bound or not.
                                  This should only be set to False for demonstration or testing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_sigma_embeddings(z)
            [1.306, 1.377]
            sage: P1._bound_for_sigma_embeddings(z,use_initial_bd_d=False)
            [1.306, 1.377]
            sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,0.5)])
            sage: P1._bound_for_sigma_embeddings(z)
            [1.273, 1.273]
            sage: P1._bound_for_sigma_embeddings(z,use_initial_bd_d=False)
            [2.545, 2.545]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_sigma_embeddings(z)
            [3.605, 3.800]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P4._bound_for_sigma_embeddings(z)
            [1.984, 2.091]

        """
        self._check_upper_half_plane_element(z)
        d = self._Dz(z, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d)
        return [upper(d[i] * y ** (-1 / 2), prec=prec) for i, y in enumerate(z.imag())]

    def _bound_for_sigma_coordinates(self, z, initial_bd_d=None, prec=16, use_initial_bd_d=True):
        """
        Bound `c` for the coordinates, with respect to the ring of integers,
        of the closest cusp to z.


        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``initial_bd_d`` -- float - an initial bound (default None) for the distance to a cusp.
                              If it is None or larger than the distance to infinity
                              then this distance is used.
        - ``prec`` -- the number of bits precision in the returned values.
        - ``use_initial_bd_d`` -- boolean (default: 'True') Use the initial bound or not.
                                  This should only be set to False for demonstration or testing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_sigma_coordinates(z)
            [1.358, 1.200]

            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_sigma_coordinates(z)
            [3.702, 1.171]

            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P4._bound_for_sigma_coordinates(z)
            [2.038, 1.177]


        """
        self._check_upper_half_plane_element(z)
        n = self.group().number_field().degree()
        d = self._Dz(z, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d)
        bounds = []
        B = self.basis_matrix_ideal().inverse()
        for i in range(n):
            bd = 0
            for j, y in enumerate(z.imag()):
                bd += B[i, j].abs() * y ** (-1 / 2)
            bounds.append(upper(bd * d[i], prec=prec))
        return bounds

    def _bound_for_rho_embeddings(
        self, z, sigma, initial_bd_d=None, use_initial_bd_d=True, prec=16
    ):
        """
        Bound for the embeddings of the numerator of the closest cusp to z.


        ``rho_i in [sigma_i x_i - delta*y_i**0.5 , sigma_i x_i - delta*y_i**0.5]``

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``sigma`` -- list of complex embeddings of algebraic integer
        - ``prec`` -- the number of bits precision in the returned values.
        - ``initial_bd_d`` -- an initial bound for the distance to nearest cusp (default None)
        - ``use_initial_bd_d`` -- boolean (default: 'True') Use the initial bound or not.
                                  This should only be set to False for demonstration or testing.


        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_rho_embeddings(z, 0)
            [(-1.306, 1.306), (-1.239, 1.239)]
            sage: b1,b2=H1.number_field().OK().gens();
            sage: P1._bound_for_rho_embeddings(z, b1)
            [(-1.925, 0.6880), (0.3790, 2.857)]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_rho_embeddings(z, 0)
            [(-3.605, 3.605), (-3.420, 3.420)]
            sage: b1,b2 = H2.number_field().OK().gens();
            sage: P2._bound_for_rho_embeddings(z, b2)
            [(-6.767, 0.4421), (-0.2571, 6.582)]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P4._bound_for_rho_embeddings(z, 0)
            [(-1.984, 1.984), (-1.882, 1.882)]
            sage: b1,b2=H4.number_field().OK().gens();
            sage: P4._bound_for_rho_embeddings(z, b1)
            [(-0.9835, 2.984), (-0.8817, 2.882)]

        """
        self._check_upper_half_plane_element(z)
        n = self.group().number_field().degree()
        d = self._Dz(z, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d)
        if not isinstance(sigma, list):
            sigma = self.group().number_field()(sigma)
            sigma = sigma.complex_embeddings()
        bounds = []
        for i in range(n):
            y = z.imag()[i]
            xs = z.real()[i] * sigma[i]
            dy = d[i] * y ** (0.5)
            b0 = xs - dy
            b1 = xs + dy
            bounds.append((b0, b1))
        res = []
        for b0, b1 in bounds:
            # We bound the lower bound differently depending on whether it is positive or negative
            b0 = lower(b0, prec=prec)
            b1 = upper(b1, prec=prec)
            res.append((b0, b1))
        return res

    def _bound_for_rho_coordinates(
        self, z, sigma, initial_bd_d=None, use_initial_bd_d=True, prec=16
    ):
        """
        Bound for the coordinates, with respect to the ring of integers,
        of the numerator of the closest cusp to z.


        ``rho_i in [sigma_i x_i - delta*y_i**0.5 , sigma_i x_i - delta*y_i**0.5]``

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``sigma`` -- complex embeddings of algebraic integer
        - ``prec`` -- the number of bits precision in the returned values.
        - ``initial_bd_d`` -- an initial bound for the distance to nearest cusp (default None)
        - ``use_initial_bd_d`` -- boolean (default: 'True') Use the initial bound or not.
                                This should only be set to False for demonstration or testing.


        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_rho_coordinates(z, 0)
            [2.571, 2.571]
            sage: b1,b2=H1.OK().gens();
            sage: P1._bound_for_rho_coordinates(z, b1)
            [5.839, 4.205]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_rho_coordinates(z, 0)
            [7.094, 7.094]
            sage: b1,b2=H2.OK().gens();
            sage: P2._bound_for_rho_coordinates(z, b2)
            [20.39, 20.39]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal = lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P4._bound_for_rho_coordinates(z, 0)
            [3.904, 3.904]
        """
        self._check_upper_half_plane_element(z)
        n = self.group().number_field().degree()
        d = self._Dz(z, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d)
        if not isinstance(sigma, list):
            sigma = self.group().number_field()(sigma)
            sigma = sigma.complex_embeddings()
        bounds = []
        factor = 1.01
        sz = [(sigma[j] * z.real()[j]).abs() for j in range(n)]
        dy = [d[j] * z.imag()[j] ** 0.5 for j in range(n)]
        for i in range(n):
            bd = 0
            for j, _y in enumerate(z.imag()):
                bd += dy[j] + sz[j] * self.basis_matrix_ideal()[i, j].abs()

            bounds.append(upper(bd * factor, prec=prec))
        return bounds

    def _candidate_integers_sigma(
        self,
        z,
        domain="polytope",
        return_polyhedron=False,
        ideal_basis=None,
        lattice_basis=None,
        do_sort=True,
        initial_bd_d=None,
        use_initial_bd_d=True,
        use_norm_bound=True,
    ):
        """
        Compute a list of candidates for the denominator of the closest cusp.

        INPUT:

        - ``z`` -- element of type UpperHalfPlaneProductelement_class
        - ``a`` -- ideal or algebraic integer (default = 1).
                   If an integer is given then the ideal is the principal ideal.
        - ``domain`` -- string: 'polytope' (default), 'boundingbox',
        - ``return_polyhedron`` -- boolean
        - ``ideal_basis`` -- list or = None,
        - ``lattice_basis`` -- list of lists corresponding to a numerical basis of the lattice
                               corresponding to the ring of integers.
        - ``sorted`` -- boolean -- True to return a list sorted by norm first and then
                                   lexicographically with respect to embeddings.
        - ``initial_bd_d`` -- positive number (default: `None`) - an initial bound for the distance
                                                                  to the nearest cusp.
        - ``use_norm_bound`` -- boolean (default: `True`) -- True if using the norm bound otherwise
        - ``use_initial_bd_d`` -- boolean (default: 'True') Use the initial bound or not.
                                  This should only be set to False for demonstration or testing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1._candidate_integers_sigma(z)
            [0, -1, 1]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1._candidate_integers_sigma(z)
            [0, -1, 1]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.7)])
            sage: P1._candidate_integers_sigma(z)
            [0, -1, 1]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P2._candidate_integers_sigma(z,use_norm_bound=False)[0:9]
            [0, -1, 1, -2, 2, -3, 3, a, -a]
            sage: P2._candidate_integers_sigma(z)
            [0, -1, 1]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: P4._candidate_integers_sigma(z)
            [0, -1, 1]
        """
        if return_polyhedron:
            if domain == "boundingbox":
                bounds = self._bound_for_sigma_coordinates(
                    z, initial_bd_d=initial_bd_d, use_initial_bd=use_initial_bd_d
                )
                B = None
            elif domain == "preimage":
                bounds = self._bound_for_sigma_embeddings(
                    z, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d
                )
                B = None
            else:
                bounds = self._bound_for_sigma_embeddings(
                    z, initial_bd_d=initial_bd_d, use_initial_bd_d=use_initial_bd_d
                )
                B = self.basis_matrix_ideal().inverse()
            return self.polytope_from_bounds(bounds, B)
        # Else we use efficient methods to find candidates.
        sigma_candidates = find_candidate_cusps(
            self,
            z,
            return_sigma_candidates=True,
            use_norm_bound=use_norm_bound,
            initial_bd_d=initial_bd_d,
            use_initial_bd_d=use_initial_bd_d,
        )
        if do_sort:

            def absort(val):
                return (sum([abs(x) ** 2 for x in val.complex_embeddings()]),) + tuple(
                    val.complex_embeddings()
                )

            sigma_candidates.sort(key=absort)
        return sigma_candidates

    def _get_lattice_and_ideal_basis(self, ideal=None):
        """
        Compute an integral basis for an ideal as well as a basis matrix for the associated lattice
         in R^n in the form of a nested list.

        INPUT:

        - ``ideal`` -- ideal (default: `None`) if none then (1) is used.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: P1._get_lattice_and_ideal_basis()
            ([[1.00000000000000, -1.61803398874989],
             [1.00000000000000, 0.618033988749895]],
             [1, 1/2*a - 1/2])
            sage: P1._get_lattice_and_ideal_basis(H1.OK().fractional_ideal(2))
            ([[2.00000000000000, -3.23606797749979], [2.00000000000000, 1.23606797749979]],
            [2, a - 1])
            sage: H2=ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: P2._get_lattice_and_ideal_basis()
             ([[1.00000000000000, -3.16227766016838], [1.00000000000000, 3.16227766016838]],
              [1, a])
            sage: a= H2.OK().fractional_ideal(H2.OK().gen(1)+1)
            sage: P2._get_lattice_and_ideal_basis(a)
            ([[9.00000000000000, -2.16227766016838], [9.00000000000000, 4.16227766016838]], [9, a + 1])
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: P4._get_lattice_and_ideal_basis()
            ([[1.00000000000000, -1.73205080756888], [1.00000000000000, 1.73205080756888]],
             [1, a])
            sage: a = H4.OK().fractional_ideal(H4.OK().gen(1)+1)
            sage: P4._get_lattice_and_ideal_basis(a)
            ([[2.00000000000000, -0.732050807568877],
               [2.00000000000000, 2.73205080756888]],
               [2, a + 1])



        """
        if not ideal:
            ideal = self._construct_ideal(1)
        ideal_basis = ideal.basis()
        lattice_basis = self.basis_matrix_ideal(ideal)
        n = len(lattice_basis[0])
        # Make lattice basis to a nested list
        # to avoid creation of FreeModule elements
        lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
        return lattice_basis, ideal_basis

    def _candidate_integers_rho(
        self,
        z,
        sigma,
        a=1,
        domain="polytope",
        return_polyhedron=False,
        ideal_basis=None,
        lattice_basis=None,
        do_sort=True,
        use_initial_bd_d=True,
    ):
        """
        Compute a list of candidates for the denominator of the closest cusp.

         INPUT:

        - ``z`` -- element of type UpperHalfPlaneProductelement_class
        - ``sigma`` -- algebraic integer
        - ``a`` -- ideal or algebraic integer (default = 1).
                  If an integer is given then the ideal is the principal ideal.
        - ``algorithm`` -- string (either 'coordinates' or 'embeddings')
        - ``use_initial_bd_d`` -- boolean (default: 'True')
                  Use the initial bound or not.
                  This should only be set False for demonstration or testing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1._candidate_integers_rho(z,1)
            [0, 1, 2]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.25)])
            sage: P1._candidate_integers_rho(z,1)
            [1]
            sage: H2 = ExtendedHilbertModularGroup(10)
            sage: P2 = ExtendedHilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P2._candidate_integers_rho(z,1)
            [0, -1, 1, -2, 2, 3, a + 1, -a + 1, 4]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P4._candidate_integers_rho(z, 1)
            [0, 1, a + 1, 2, -a + 1]
        """
        if return_polyhedron:
            if domain == "boundingbox":
                bounds = self._bound_for_rho_coordinates(z, sigma)
                B = None
            elif domain == "preimage":
                bounds = self._bound_for_rho_embeddings(z, sigma)
                B = None
            else:
                bounds = self._bound_for_rho_embeddings(z, sigma)
                B = self.basis_matrix_ideal().inverse()
            return self.polytope_from_bounds(bounds, B)
        if not lattice_basis or not ideal_basis:
            lattice_basis, ideal_basis = self._get_lattice_and_ideal_basis()
        dist = None
        if use_initial_bd_d:
            candidate_cusp = self.get_heuristic_closest_cusp(z)
            if candidate_cusp:
                dist = distance_to_cusp_eg(self, candidate_cusp[0], candidate_cusp[1], z)
        rho_coordinate_bounds = self._bound_for_rho_coordinates(
            z, sigma, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d
        )
        rho_coordinate_bounds = [(-b, b) for b in rho_coordinate_bounds]
        rho_embedding_bounds = self._bound_for_rho_embeddings(
            z, sigma, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d
        )
        rho_candidates_coordinates = lattice_elements_in_box(
            lattice_basis, rho_embedding_bounds, rho_coordinate_bounds
        )
        rho_candidates = coordinates_to_ideal_elements(rho_candidates_coordinates, ideal_basis)
        if do_sort:

            def absort(val):
                return (sum([abs(x) ** 2 for x in val.complex_embeddings()]),) + tuple(
                    val.complex_embeddings()
                )

            rho_candidates.sort(key=absort)
        return rho_candidates

    def _candidate_closest_cusps(
        self, z, use_lll=True, use_norm_bound=True, use_initial_bd_d=True, as_cusps=False
    ):  # Doing
        r"""
        Find candidates for the closest cusp.

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``use_lll`` -- boolean (default: `True`)
                                Use the LLL method to find a preliminary bounds
        - ``use_norm_bound`` -- boolean (default: `True`)
                                Use the norm bound together with the embedding
                                bounds
        - ``use_initial_bd_d`` -- boolean (default: `False`) Use initial bound

        - ``as_cusps`` -- boolean - (default: `False`), set to True to return a
                                            list of cusps instead of tuples.


        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import *
            sage: from hilbert_modgroup.all import *
            sage: H1 = ExtendedHilbertModularGroup(5)
            sage: P1 = ExtendedHilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0, 1),CC(0, 1)])
            sage: P1._candidate_closest_cusps(z)
            [(1, 0), (-1, -1), (0, 1), (1, -1)]
            sage: K = QuadraticField(3)
            sage: lattice_ideal = K.different()
            sage: H4 = ExtendedHilbertModularGroup(K, lattice_ideal)
            sage: P4 = ExtendedHilbertPullback(H4)
            sage: z = UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P4._candidate_closest_cusps(z)
            [(1, 0), (-2, -1), (-a - 1, -1), (-1, -1), (a - 1, -1), (0, 1)]
        """
        cusp_candidates = find_candidate_cusps(
            self,
            z,
            use_lll=use_lll,
            use_norm_bound=use_norm_bound,
            return_sigma_candidates=False,
            use_initial_bd_d=use_initial_bd_d,
        )
        if as_cusps:
            # Convert to cusps
            cusps = []
            lattice_ideal = self.group().lattice_ideal()
            for rho, sigma in cusp_candidates:
                c = NFCusp_wrt_lattice_ideal(lattice_ideal, rho, sigma)
                if c not in cusps:
                    cusps.append(c)
            return cusps
        else:
            return cusp_candidates

    def level_reduction_matrix(self, M):
        r"""
        Return a coset representative ``A`` such that ``A * M`` belongs to the group.

        Given a matrix ``M``, this method searches through the right coset
        representatives of the group and returns a matrix ``A`` such that
        ``A * M`` is an element of the extended Hilbert modular group.

        INPUT:

        - ``M`` -- a matrix

        OUTPUT:

        A matrix ``A`` from the coset representatives such that ``A * M``
        is in the group.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
            sage: K.<a> = QuadraticField(5)
            sage: level_ideal = K.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K, level_ideal=level_ideal)
            sage: P = ExtendedHilbertPullback(H)
            sage: M = H(1).matrix()
            sage: A = P.level_reduction_matrix(M)
            sage: A * M in H
            True

        """
        coset_matrices = self.group().coset_matrices()
        Mat = None
        for A in coset_matrices:
            if A * M in self.group():
                Mat = A
                break
        if Mat is None:
            raise ValueError("No coset representative B satisfies B*M ∈ G.")
        return Mat

    def reduce(self, z, return_map=False):
        r"""
        Reduce ``z`` to a point in the fundamental domain.

        INPUT:

        - ``z`` -- point in the upper half-plane
        - ``return_map`` -- boolean (default False)
          Set to ``True`` to return the map which performed the reduction.

           EXAMPLES::
              sage: from hilbert_modgroup.all import *
              sage: from hilbert_modgroup.extended.all import *
              sage: K5.<a> = QuadraticField(5)
              sage: H5 = ExtendedHilbertModularGroup(K5)
              sage: P5 = ExtendedHilbertPullback(H5)
              sage: H5L = HilbertModularGroup(K5)
              sage: P5L = HilbertPullback(H5L)
              sage: z = UpperHalfPlaneProductElement([0.025+I/2,1.0233+I])
              sage: P5.reduce(z)
              [-0.656135336505516 + 0.762026955112429*I, 0.309364137820355 + 1.27887213029370*I]
              sage: P5L.reduce(z)
              [-0.656135336505516 + 0.762026955112429*I, 0.309364137820355 + 1.27887213029370*I]
              sage: level_ideal = K5.fractional_ideal(2)
              sage: H5 = ExtendedHilbertModularGroup(K5, level_ideal = level_ideal)
              sage: P5 = ExtendedHilbertPullback(H5)
              sage: P5.reduce(z)
              [-1.03996773704886 + 0.345717054166798*I, 0.257798072768767 + 0.239024655533354*I]
              sage: lattice_ideal = K5.different()
              sage: H5 = ExtendedHilbertModularGroup(K5, lattice_ideal = lattice_ideal, level_ideal = level_ideal)
              sage: P5 = ExtendedHilbertPullback(H5)
              sage: P5.reduce(z)
              [0.472213595499958 + 0.500000000000000*I, 0.576086404500042 + 1.00000000000000*I]
              sage: K3.<a> = QuadraticField(3)
              sage: H3 = ExtendedHilbertModularGroup(K3)
              sage: H3L = HilbertModularGroup(K3)
              sage: P3 = ExtendedHilbertPullback(H3)
              sage: P3L = HilbertPullback(H3L)
              sage: P3.reduce(z)
              [-0.758778906564500 + 0.534561979912464*I, 0.866520262514823 + 1.82305340081506*I]
              sage: P3L.reduce(z)
              [-0.0997506234413965 + 1.99501246882793*I, -0.499867403002826 + 0.488485686507208*I]
              sage: H3 = ExtendedHilbertModularGroup(K3, tp_units = False)
              sage: P3 = ExtendedHilbertPullback(H3)
              sage: P3.reduce(z)
              [-0.0997506234413965 + 1.99501246882793*I, -0.499867403002826 + 0.488485686507208*I]


        """
        K = self.number_field()
        lattice_ideal = self.group().lattice_ideal()
        level_ideal = K.fractional_ideal(1)
        tp_units = self.group().tp_units()
        H = ExtendedHilbertModularGroup(
            K, lattice_ideal=lattice_ideal, level_ideal=level_ideal, tp_units=tp_units
        )
        P = ExtendedHilbertPullback(H)
        if z.norm() == 0:
            raise ValueError("Can not reduce point at the boundary of one of the half-planes.")
        c = P.find_closest_cusp(z, return_multiple=False, as_cusp=True)
        c_rep, Umu = P.group().cusp_representative(c, return_map=True)
        A = P._group.cusp_normalizing_map(c_rep)
        w = z.apply(A.inverse() * Umu)
        w, B = P.reduce_in_cuspidal_region(w, c_rep, return_map=True)
        w = w.apply(A)
        M = A * B * A.inverse() * Umu
        Mat = self.level_reduction_matrix(M)
        w = w.apply(Mat)
        if return_map:
            return w, Mat * M
        return w
