r"""
A utility class for managing the reduction algorithm for Hilbert modular groups

AUTHORS:

- Fredrik Stromberg (2021)

NOTE: I know it is often a bad idea to write utility classes but I decided to do it anyway at least for the moment.
"""
from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup_class
from sage.categories.sets_cat import cartesian_product
from sage.functions.other import floor
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modular.cusps_nf import NFCusp
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from sage.rings.number_field.number_field_ideal import is_NumberFieldIdeal
from sage.rings.real_double import RDF
from sage.rings.real_mpfr import RealField
from sage.structure.sage_object import SageObject

from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement__class

from hilbert_modgroup.utils import upper, lower

from hilbert_modgroup.pullback_cython import lattice_elements_in_box, \
    coordinates_to_ideal_elements


class HilbertPullback(SageObject):
    r"""
    Utility class for pullback / reduction algorithms for Hilbert modular groups.

    """
    def __init__(self,G):
        r"""
        Init self.

        INPUT::

        - ``G`` - A Hilbert modular group

        EXAMPLES:

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H = HilbertModularGroup(5)
            sage: HilbertPullback(H)
            Pullback class for Hilbert Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?

        """
        if not isinstance(G,HilbertModularGroup_class):
            raise ValueError("Need a Hilbert modular group")
        self._group = G

    def __str__(self):
        r"""
        Return string representation of self.

         EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H = HilbertModularGroup(5)
            sage: P = HilbertPullback(H)
            sage: str(P)
            'Pullback class for Hilbert Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?'

        """
        return f"Pullback class for {self._group}"

    def __repr__(self):
        r"""
        Return string representation of self.

         EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H = HilbertModularGroup(5)
            sage: P = HilbertPullback(H)
            sage: P
            Pullback class for Hilbert Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?

        """
        return str(self)

    def group(self):
        """
        Return the group of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H = HilbertModularGroup(5)
            sage: P = HilbertPullback(H)
            sage: P.group()
            Hilbert Modular Group PSL(2) over Maximal Order in Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?

        """
        return self._group

    def _check_upper_half_plane_element(self, z):
        r"""
        Check if z is an element of type UpperHalfPlaneProductElement__class of the correct degree.

        INPUT:
        - `z` - potential element of a product of upper half planes.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertPullback, HilbertModularGroup, UpperHalfPlaneProductElement, ComplexPlaneProductElement
            sage: H = HilbertModularGroup(5)
            sage: P = HilbertPullback(H)
            sage: P._check_upper_half_plane_element([1,1])
            Traceback (most recent call last):
            ...
            ValueError: Need an element of type: UpperHalfPlaneProductElement__class of degree 2
            sage: z = UpperHalfPlaneProductElement([1+I,1+I])
            sage: P._check_upper_half_plane_element(z)
            True
            sage: z = UpperHalfPlaneProductElement([1+I,1+I,1+I])
            sage: P._check_upper_half_plane_element(z)
            Traceback (most recent call last):
            ...
            ValueError: Need an element of type: UpperHalfPlaneProductElement__class of degree 2
            sage: z = ComplexPlaneProductElement([1+I,1+I])
            sage: P._check_upper_half_plane_element(z)
            Traceback (most recent call last):
            ...
            ValueError: Need an element of type: UpperHalfPlaneProductElement__class of degree 2


        """
        if not isinstance(z, UpperHalfPlaneProductElement__class) or z.degree() != self.group().base_ring().degree():
            raise ValueError("Need an element of type: UpperHalfPlaneProductElement__class of degree {0}".format(
                self.group().base_ring().degree()))
        return True

    @cached_method
    def basis_matrix_logarithmic_unit_lattice(self, prec=53):
        """
        Return the Basis matrix for the logarithmic unit lattice / +-1

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1.basis_matrix_logarithmic_unit_lattice()
            [ 0.481211825059603]
            [-0.481211825059603]
            sage: H2=HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: P2.basis_matrix_logarithmic_unit_lattice()
            [ 1.81844645923207]
            [-1.81844645923207]

            sage: var('x')
            x
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3.basis_matrix_logarithmic_unit_lattice()
            [ 1.78943386474421 -4.27317838834468]
            [-3.58349750383703  1.78711898997458]
            [ 1.79406363909282  2.48605939837008]

        """
        n = self.group().base_ring().degree()
        entries = [[x.abs().log() for x in u.complex_embeddings(prec)] for u in self.fundamental_units()]
        return matrix(RealField(prec), n-1, n, entries).transpose()

    @cached_method()
    def fundamental_units(self):
        r"""
        Return fundamental the units for the group of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1.fundamental_units()
            [1/2*a - 1/2]
            sage: x = ZZ['x'].gen()
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3.fundamental_units()
            [-a, -a - 6]

        """
        return self.group().base_ring().number_field().unit_group().fundamental_units()

    def Y(self,z, return_error_estimate=False):
        r"""
        Compute the coordinate of y=Im(z)/N(y)^(1/n) with respect to the logarithmic unit lattice.

        INPUT::

        - ``z`` -- element of type UpperHalfPlaneProductElement
        - ``return_error_estimate`` -- boolean (default=False) set to True to return a tuple including an estimate of the error.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: u0,u1=H1.base_ring().number_field().unit_group().gens()
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.Y(z)
            (0.000000000000000)
            sage: P1.Y(u1*z)
            Traceback (most recent call last):
            ...
            ValueError: Point [-1.61803398874989 - 1.61803398874989*I, 0.618033988749895 + 0.618033988749895*I] not in upper half-plane!
            sage: P1.Y(u1**2*z)
            (2.00000000000000)
            sage: P1.Y(u1**-2*z)
            (-2.00000000000000)
            sage: z=UpperHalfPlaneProductElement([CC(1,3),CC(0,2)])
            sage: P1.Y(z)
            (0.421295869088362)
            sage: P1.Y(u1**2*z)
            (2.42129586908836)

            # The precision gets worse when the imaginary parts have large differences.
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,2),CC(0,1)])
            sage: P3.Y(z,return_error_estimate=True)
            ((-0.128907673154173, 0.0000882959672881647), 2.220446049250313e-16)
            sage: z=UpperHalfPlaneProductElement([CC(0,100),CC(0,200),CC(0,100)])
            sage: P3.Y(z,return_error_estimate=True)
            ((-0.128907673154173, 0.0000882959672880335), 6.38378239159465e-16)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,2),CC(0,100)])
            sage: P3.Y(z)
            (0.638975592292095, 0.680877344221728)
            sage: P3.Y(z,return_error_estimate=True)
            ((0.638975592292095, 0.680877344221728), 1.1546319456101628e-14)
            sage: CF=ComplexField(106)
            sage: z=UpperHalfPlaneProductElement([CF(0,1),CF(0,2),CF(0,100)])
            sage: P3.Y(z,return_error_estimate=True)
            ((0.6389755922920966140584010614944, 0.6808773442217317886140786738603),
             2.0214560696288428e-30)


        """
        self._check_upper_half_plane_element(z)
        normalized_imag = z/z.imag_norm()**(1/self.group().base_ring().degree())
        log_vector = matrix(vector(normalized_imag.imag_log())).transpose()
        B = self.basis_matrix_logarithmic_unit_lattice(prec=z.base_ring().prec())
        coordinate_vector = B.solve_right(log_vector, check=False)
        if return_error_estimate:
            return vector(coordinate_vector),(B*coordinate_vector-log_vector).norm(Infinity)
        return vector(coordinate_vector)

    def reduce_by_units(self,z, return_map=False):
        r"""
        Reduce the point z with respect to action of units to a point where Im(z) belongs
        to a fundamental domain for the logarithmic unit lattice.

        INPUT::

        - ``z`` -- element of type UpperHalfPlaneProductElement
        - ``return_map`` -- boolean (default=False) set to True to return the map which does the reduction.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: P1.reduce_by_units(z)
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(0,1000)]);z
            [1.00000000000000 + 1.00000000000000*I, 1000.00000000000*I]
            sage: P1.reduce_by_units(z,return_map=True)
            (
            [46.9787137637478 + 46.9787137637478*I, 21.2862362522082*I],
            [-3/2*a + 7/2            0]
            [           0  3/2*a + 7/2]
            )
            sage: w,A=_
            sage: z.apply(A) == w
            True

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,2),CC(0,1)]); z
            [1.00000000000000*I, 2.00000000000000*I, 1.00000000000000*I]
            sage: P3.reduce_by_units(z)
            [1.00000000000000*I, 2.00000000000000*I, 1.00000000000000*I]
            sage: z=UpperHalfPlaneProductElement([CC(0,0.001),CC(0,2),CC(0,1000)])
            sage: P3.reduce_by_units(z,return_map=True)
            (
            [5.14796503476537*I, 0.0560735644527675*I, 6.92845248925837*I],
            <BLANKLINE>
            [-a^2 + 6*a          0]
            [        0     -a - 6]
            )

            sage: z=UpperHalfPlaneProductElement([CC(0,0.0001),CC(0,0.1),CC(0,100)])
            sage: P3.reduce_by_units(z,return_map=True)
            (
            [0.0143665696311556*I, 3.63341121163000*I, 0.0191572146738068*I],
            <BLANKLINE>
            [    a - 6         0]
            [        0 a^2 + 6*a]
            )

        """
        K = self.group().base_ring().number_field()
        units = K.unit_group().gens()[1:]  # Only include the units != -1
        # To avoid overflow it is more efficient to apply the map, e.g. compute (z*u**-k)/u**k instead of z*u**-(2k)
        floors = [-floor(y/2+1/2) for y in self.Y(z)]
        reducing_map = prod([self.group().E(u ** y) for u,y in zip(units,floors)])
        reduced_point = z.apply(reducing_map)
        if return_map:
            return reduced_point, reducing_map
        return reduced_point

    def is_reduced_by_units(self,z):
        r"""
        Return True if z is reduced with respect to the logarithmic unit lattice, otherwise return False.

        INPUT::

        - ``z`` -- element of the type UpperHalfPlaneProductElement_class

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: P1.is_reduced_by_units(z)
            True
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(0,1000)]);z
            [1.00000000000000 + 1.00000000000000*I, 1000.00000000000*I]
            sage: P1.is_reduced_by_units(z)
            False
            sage: P1.reduce_by_units(z,return_map=True)
            (
            [46.9787137637478 + 46.9787137637478*I, 21.2862362522082*I],
            <BLANKLINE>
            [-3/2*a + 7/2            0]
            [           0  3/2*a + 7/2]
            )
            sage: w,A=_
            sage: P1.is_reduced_by_units(w)
            True

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,2),CC(0,1)])
            sage: P3.is_reduced_by_units(z)
            True

            sage: z=UpperHalfPlaneProductElement([CC(0,0.0001),CC(0,0.1),CC(0,100)])
            sage: P3.is_reduced_by_units(z)
            False
            sage: w,A=P3.reduce_by_units(z,return_map=True)
            sage: P3.is_reduced_by_units(w)
            True

        """
        return all(-1 <= y < 1 for y in self.Y(z))

    @cached_method
    def basis_matrix_ideal(self, a=None, prec=53):
        r"""
        Return the Basis matrix corresponding to an integer basis of an ideal a.

        INPUT::

        - ``a`` -- ideal or number field element.
        - ``prec`` -- integer (default=53)

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1.basis_matrix_ideal()
            [ 1.00000000000000 -1.61803398874989]
            [ 1.00000000000000 0.618033988749895]
            sage: P1.basis_matrix_ideal(2)
            [ 2.00000000000000 -3.23606797749979]
            [ 2.00000000000000  1.23606797749979]

            sage: H2=HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: P2.basis_matrix_ideal()
            [ 1.00000000000000 -3.16227766016838]
            [ 1.00000000000000  3.16227766016838]
            sage: a=H2.base_ring().gen(1)
            sage: P2.basis_matrix_ideal(a+1)
            [ 9.00000000000000 -2.16227766016838]
            [ 9.00000000000000  4.16227766016838]

            sage: var('x')
            x
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3.basis_matrix_ideal()
            [   1.00000000000000   -5.98606258583498    2.28229423189948]
            [   1.00000000000000 -0.0277783731902446   -7.67566891172438]
            [   1.00000000000000    6.01384095902523    6.39337467982490]
            sage: c = P3.group().base_ring().class_group()[1].ideal()
            sage: P3.basis_matrix_ideal(c)
            [ 2.00000000000000 -4.98606258583498  3.28229423189948]
            [ 2.00000000000000 0.972221626809755 -6.67566891172438]
            [ 2.00000000000000  7.01384095902523  7.39337467982490]

        """
        ideala = self._construct_ideal(a)
        entries = [[x for x in beta.complex_embeddings(prec)] for beta in ideala.integral_basis()]
        n = self.group().base_ring().degree()
        return matrix(RealField(prec), n, n, entries).transpose()

    @cached_method
    def basis_matrix_ideal__norm(self, a=None, prec=53, row=None):
        r"""
        Return the Basis matrix corresponding to an integer basis of an ideal a.

        INPUT::

        - ``a`` -- ideal or number field element.
        - ``prec`` -- integer (default=53)

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1.basis_matrix_ideal__norm() # abs tol 1e-10
            2.61803398874989
            sage: P1.basis_matrix_ideal__norm(2) # abs tol 1e-10
            5.23606797749979
            sage: H2=HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: P2.basis_matrix_ideal__norm() # abs tol 1e-10
            4.16227766016838
            sage: a=H2.base_ring().gen(1)
            sage: P2.basis_matrix_ideal__norm(a+1) # abs tol 1e-10
            13.16227766016838
            sage: x = ZZ['x'].gen()
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3.basis_matrix_ideal__norm() # abs tol 1e-10
            13.40721563885013
            sage: c = P3.group().base_ring().class_group()[1].ideal()
            sage: P3.basis_matrix_ideal__norm(c) # abs tol 1e-10
            16.407215638850133

        """
        B = self.basis_matrix_ideal(a,prec=prec)
        if row is None:
            return B.norm(Infinity)
        elif 0 <= row < B.nrows():
            return sum([abs(x) for x in B.row(row)])
        else:
            raise ValueError(f"Can not find row:{row}")

    def _construct_ideal(self, a):
        r"""
        Construct an ideal of the number field associated with self from a

        EXAMPLES::

        sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
        sage: P1=HilbertPullback(HilbertModularGroup(5))
        sage: P1._construct_ideal(1)
        Fractional ideal (1)
        sage: P1._construct_ideal(P1.number_field().ideal(1))
        Fractional ideal (1)


        """
        if a is None:
            ideala = self.group().base_ring().ideal(1)
        elif is_NumberFieldIdeal(a):
            ideala = a
        elif a in self.group().base_ring():
            ideala = self.group().base_ring().ideal(a)
        else:
            raise ValueError("Could not construct a number field ideal from a={0}".format(a))
        return ideala

    def X(self, z, a=None):
        r"""
        Coordinate of z with respect to the integral basis of an ideal.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.X(z)
            (1.00000000000000, 0.000000000000000)
            sage: b1,b2 = H1.base_ring().ideal(1).integral_basis(); b1,b2
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

    def reduce_by_translations(self,z, a=None, return_map=False):
        r"""
        Reduce the point z with respect to the cuspidal region in a neighbourhood of the cusp .

        INPUT::

        - ``z`` -- point in the upper half-plane.
        - ``a`` -- ideal cusp (default=None) if None use the ring of integers.
        - ``return_map`` -- boolean (default=False), if set to True also return the map that makes the reduction.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.reduce_by_translations(z,return_map=True)
            (
                                                      [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )


            sage: z=UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: P1.reduce_by_translations(z,return_map=True) # abs tol 1e-10
            (
            [-0.236067977499790 + 1.00000000000000*I, 0.236067977499790 + 1.00000000000000*I],
            <BLANKLINE>
            [    1 a - 1]
            [    0     1]
            )
            sage: b1,b2=H1.base_ring().gens(); b1+b2
            3/2*a + 1/2
            sage: P1.reduce_by_translations(z,b1+b2,return_map=True) # abs tol 1e-10
            (
            [4.76393202250021 + 1.00000000000000*I, 5.23606797749979 + 1.00000000000000*I],
            <BLANKLINE>
            [    1 a + 4]
            [    0     1]
            )

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P3.reduce_by_translations(z,return_map=True)
            Traceback (most recent call last):
            ...
            ValueError: Need an element of type: UpperHalfPlaneProductElement__class of degree 3
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1),CC(1,1)])
            sage: P3.reduce_by_translations(z,return_map=True)
            (
                                                                          [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(2,1),CC(10,1)])
            sage: P3.reduce_by_translations(z,return_map=True)
            (
            [2.98606258583498 + 1.00000000000000*I, -1.97222162680976 + 1.00000000000000*I, -0.0138409590252291 + 1.00000000000000*I],
            <BLANKLINE>
            [     1 -a - 4]
            [     0      1]
            )
            sage: w,A=_
            sage: vector([floor(x+0.5) for x in P3.X(z)]) + A[0,1].vector() == 0
            True
        """
        X = self.X(z, a)
        ideala = self._construct_ideal(a)
        basis = ideala.integral_basis()
        correction = sum([b*floor(X[i]+0.5) for i, b in enumerate(basis)])
        reduced_point = z - correction
        if return_map:
            return reduced_point, self.group().T(-correction)
        return reduced_point

    def is_reduced_by_translations(self, z, a=None, prec=53):
        r"""
        Check if the given point is reduced in the cuspidal region.

        EXAMPLES::


            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.is_reduced_by_translations(z)
            False
            sage: P1.reduce_by_translations(z)
            [1.00000000000000*I, 1.00000000000000*I]
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z))
            True
            sage: z=UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: P1.is_reduced_by_translations(z)
            False
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z))
            True
            sage: b1,b2=H1.base_ring().gens(); b1+b2
            3/2*a + 1/2
            sage: P1.is_reduced_by_translations(z,b1+b2)
            False
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z,b1+b2),b1+b2)
            True
            sage: P1.is_reduced_by_translations(P1.reduce_by_translations(z,b1+b2))
            False
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(2,1),CC(10,1)])
            sage: P3.is_reduced_by_translations(z)
            False
            sage: P3.is_reduced_by_translations(P3.reduce_by_translations(z))
            True

        """
        X = self.X(z, a)
        return all(-1/2 <= x < 1/2 for x in X)

    def reduce(self, z, return_map=False):
        r"""
        Reduce ``z`` to a point in the fundamental domain.

        INPUT::

        - ``z`` -- point in the upper half-plane
        - ``return_map`` -- boolean (default False) set to ``True`` to return the map which performed the reduction.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1=HilbertModularGroup(5)
            sage: P1=HilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([1+I,1+I])
            sage: P1.reduce(z)
            [1.00000000000000*I, 1.00000000000000*I]
            sage: z = UpperHalfPlaneProductElement([0.25+I/2,1+I])
            sage: P1.reduce(z) # abs tol 1e-10
            [0.694427190999916 + 0.611145618000168*I, -0.309016994374947 + 1.30901699437495*I]
            sage: P1.reduce(z, return_map=True)[1]
            [ 1/2*a - 1/2            0]
            [-1/2*a + 1/2  1/2*a + 1/2]

        """
        c = self.find_closest_cusp(z)
        if return_map:
            A = self._group.cusp_normalizing_map(c)
            w, B = self.reduce_in_cuspidal_region(z, c, return_map=return_map)
            return w,A*B*A.inverse()
        return self.reduce_in_cuspidal_region(z, c)


    def reduce_in_cuspidal_region(self,z, cusp=None, check=True, return_map=False):
        r"""
        Reduce the point z with respect to the cuspidal region in a neighbourhood of the cusp .

        INPUT::

        - ``z`` -- point in the upper half-plane.
        - ``cusp`` -- cusp (default=None) if none using the cusp at infinity, otherwise reduce with respect to this cusp.
        - ``return_map`` -- boolean (default=False) return the map A such that AN^-1z is reduced where N is the cusp
                            normalizing map for the cusp.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.reduce_in_cuspidal_region(z)
            [1.00000000000000*I, 1.00000000000000*I]
            sage: P1.reduce_in_cuspidal_region(z,return_map=True)
            (
                                                      [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )
            sage: P1.reduce_in_cuspidal_region(z,cusp=NFCusp(P1.group().base_ring().number_field(),0,1),return_map=True)
            (
            [-0.500000000000000 + 0.500000000000000*I, -0.500000000000000 + 0.500000000000000*I],
            <BLANKLINE>
            [1 0]
            [0 1]
            )

            # Check that if we apply a cus-normalizing map to a point then the reduction with respect to that cusp is the same
            # as the reduction of the original point wih respect to infinity
            sage: c=NFCusp(P1.group().base_ring().number_field(),1,2)
            sage: N=H1.cusp_normalizing_map(c)
            sage: w=z.apply(N); w
            [0.384615384615385 + 0.0769230769230769*I, 0.384615384615385 + 0.0769230769230769*I]
            sage: P1.reduce_in_cuspidal_region(w,cusp=c,return_map=True)
            (
                                                      [ 1 -1]
            [1.00000000000000*I, 1.00000000000000*I], [ 0  1]
            )

            sage: z=UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: P1.reduce_in_cuspidal_region(z)
            [-0.236067977499790 + 1.00000000000000*I, 0.236067977499790 + 1.00000000000000*I]
            sage: P1.reduce_in_cuspidal_region(z,return_map=True)
            (
            [-0.236067977499790 + 1.00000000000000*I, 0.236067977499790 + 1.00000000000000*I],
            <BLANKLINE>
            [     1 a - 1]
            [     0      1]
            )
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(2,1),CC(10,1)])
            sage: P3.reduce_in_cuspidal_region(z,return_map=True)
            (
            [2.98606258583498 + 1.00000000000000*I, -1.97222162680976 + 1.00000000000000*I, -0.0138409590252291 + 1.00000000000000*I],
            <BLANKLINE>
            [     1 -a - 4]
            [     0      1]
            )


        """
        self._check_upper_half_plane_element(z)
        cusp_infinity = NFCusp(self.group().base_ring().number_field(), Infinity)
        if cusp is None:
            cusp = cusp_infinity
        # First normalize with respect to the cusp
        if cusp != cusp_infinity:
            N = self.group().cusp_normalizing_map(cusp, inverse=True)
            z = z.apply(N)
        ideala = cusp.ideal()**-2
        # Then reduce with respect to the units, followed by reduction by translation with respect to the ideal a**-2
        if return_map:
            w, A = self.reduce_by_units(z, return_map=return_map)
            w, B = self.reduce_by_translations(w, ideala, return_map=return_map)
            return w, B*A
        else:
            w = self.reduce_by_units(z, return_map=return_map)
            w = self.reduce_by_translations(w, ideala, return_map=return_map)
        return w

    def number_field(self):
        """
        Return number field of the Hilbert modular group of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3.number_field()
            Number Field in a with defining polynomial x^3 - 36*x - 1

        """
        return self.group().base_ring().number_field()

    def find_closest_cusp(self, z, return_multiple=False):
        r"""
        Find the closest cusp (rho:sigma) to z

        INPUT::

        - `z` -- point in the upper half-plane
        - `return_multiple` -- boolean: default False - set to True to return all cusps with the same minimal distance ,
                                        otherwise just return one closest cusp.


        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1.find_closest_cusp(z)
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
            sage: P1.find_closest_cusp(z,return_multiple=True)
            [Cusp Infinity of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?,
             Cusp [-1: -1] of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?]
            sage: H2 = HilbertModularGroup(10)
            sage: K2 = H2.base_ring().number_field()
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(0,1.0),CC(0,1.0)])
            sage: P2.find_closest_cusp(z)
            Cusp Infinity of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380?
            sage: P2.find_closest_cusp(z,return_multiple=True)
            [Cusp Infinity of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380?,
            Cusp [0: 1] of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380?]
            sage: z=UpperHalfPlaneProductElement([CC(2.58,0.5),CC(0.5,0.5)])
            sage: P2.find_closest_cusp(z)
            Cusp [3*a - 10: a - 4] of Number Field in a with defining polynomial x^2 - 10 with a = 3.162277660168380?

        """
        candidates = self._candidate_closest_cusps(z)
        min_distance = self.distance_to_cusp((1,0),z)
        cusp = NFCusp(self.number_field(), 1, 0)
        if return_multiple:
            closest_cusp = []
        else:
            closest_cusp = cusp
        for cusp in candidates:
            if not isinstance(cusp,NFCusp):
                cusp = NFCusp(self.number_field(),cusp[0],cusp[1])
            d = self.distance_to_cusp(cusp,z)
            if d < min_distance:
                min_distance = d
                if return_multiple:
                    closest_cusp = [cusp]
                else:
                    closest_cusp = cusp
            elif d == min_distance and return_multiple:
                closest_cusp.append(cusp)
            if d < self._bound_for_closest_cusp():
                break
        return closest_cusp

    def distance_to_cusp(self, cusp, z):
        """
        Give the distance from the point z to the cusp: N(a)^-1 N(Im A^-1z)^(-1/2)

        INPUT::
        - `cusp` -- NF cusp or tuple of integral elements in the number field.
        - `z`

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, ComplexPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=ComplexPlaneProductElement([CC(0,1),CC(0,1)]);
            sage: P1.distance_to_cusp(H1.cusps()[0],z)
            1.00000000000000
            sage: P1.distance_to_cusp(H1.cusps()[0],z*[2,2])
            0.500000000000000

            sage: H2=HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=ComplexPlaneProductElement([CC(0,1),CC(0,1)]);
            sage: P2.distance_to_cusp(H2.cusps()[0],z)
            1.00000000000000
            sage: P2.distance_to_cusp(H2.cusps()[1],z)
            6.40312423743285
            sage: c=NFCusp(H2.base_ring().number_field(),3,H2.base_ring().number_field().gen()+1)
            sage: P2.distance_to_cusp(c,z)
            6.32455532033675
            sage: z1=ComplexPlaneProductElement([CC(-1.38,0.1),CC(0.72,0.1)]); z1
            [-1.38000000000000 + 0.100000000000000*I, 0.720000000000000 + 0.100000000000000*I]
            sage: P2.distance_to_cusp(H2.cusps()[0],z1)
            10.0000000000000
            sage: c=NFCusp(H2.base_ring().number_field(),3,H2.base_ring().number_field().gen()+1)
            sage: P2.distance_to_cusp(c,z1)
            0.300834689631305
            sage: P2.distance_to_cusp(H2.cusps()[1],z1)
            4.12502547927946

        """
        if not isinstance(cusp,NFCusp) and not isinstance(cusp,tuple):
            raise TypeError(f"Can not convert {cusp} to a NFCusp")
        if not isinstance(cusp, NFCusp):
            cusp = NFCusp(self.group().base_ring().number_field(),cusp[0],cusp[1])

        w = self._group.apply_cusp_normalizing_map(cusp,z,inverse=True)
        return (w.imag_norm())**(-1/2)*cusp.ideal().norm()**-1

    def polytope_from_bounds(self, bounds, B=None):
        r"""
        Return the polytope defined by a set of bounds, either after applying the map 'B', or not.

        Get a list of integer points in an ideal with embeddings in RR^n within a given bound.
        Reference: Algorithm 10

        INPUT::

        - ``a`` -- ideal or algebraic integer.
        - ``bounds`` -- list of bounds for the coordinates or scalar which then gives a cube with the same bound.
                        can be of the form [b1,b2,...,bn] or [(a1,b1),(a2,b2)....(an,bn)]
                        in the first case the bounds are interpreted as (-b1,b1),...,(-bn,bn)
        - ``return_polyhedron`` -- bolean (default False) set to True to return a polyhedron of the corresponding domain
        - ``preimage`` -- boolean (default False) set to True to return the polyhedron of the pre-image (only used when return_polyhedron=True)

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: p=P1.polytope_from_bounds([(-1,1),(-1,1)])
            sage: p.vertices()
            (A vertex at (-1.0, -1.0),
             A vertex at (-1.0, 1.0),
             A vertex at (1.0, -1.0),
             A vertex at (1.0, 1.0))
            sage: p=P1.polytope_from_bounds([1,1])
            sage: p.vertices()
            (A vertex at (-1.0, -1.0),
             A vertex at (-1.0, 1.0),
             A vertex at (1.0, -1.0),
             A vertex at (1.0, 1.0))

        """
        from sage.all import Polyhedron, vector
        if not isinstance(bounds,(list,tuple)):
            raise ValueError("Need a list of bounds!")
        if not isinstance(bounds[0],tuple):
            bounds = [(-b,b) for b in bounds]
        # Hypercube we want embeddings in
        vertices = cartesian_product([[RDF(a),RDF(b)] for a,b in bounds]).list()
        p1 = Polyhedron(vertices,base_ring=RDF)
        if not B:
            return p1
        # Hypercube containing the integral points of the coordinates wrt the integral basis
        vertices = [vector([y for y in B * vector(x)]) for x in p1.vertices()]
        # Try to make a polyhedron of the mapped vertices.
        return Polyhedron(vertices,base_ring=RDF)

    @cached_method
    def max_ideal_norm(self):
        r"""
        Compute the maximum of the norms of all ideal representatives in the ideal class group.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1.max_ideal_norm()
            1
            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: P2.max_ideal_norm()
            3
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3.max_ideal_norm()
            16
        """
        return max([x.ideal().norm() for x in self.group().base_ring().class_group()])

    @cached_method
    def _exp_matrix_BLambda_norm(self):
        r"""
        Compute e^{|| B_{\Lambda} ||_{\infty}/2 }

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1._exp_matrix_BLambda_norm()
            1.272019649514069

            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: P2._exp_matrix_BLambda_norm()
            2.4823935345082533

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3._exp_matrix_BLambda_norm()
            20.72428345635302

        """
        return self.basis_matrix_logarithmic_unit_lattice().norm(Infinity).exp().sqrt()

    @cached_method()
    def D(self):
        r"""
        Return the bound ``D`` for this Hilbert modular group.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: P1 = HilbertPullback(HilbertModularGroup(5))
            sage: P1.D() # abs tol 1e-10
            1.272019649514069
            sage: P2 = HilbertPullback(HilbertModularGroup(10))
            sage: P2.D() # abs tol 1e-10
            4.299631726148779
            sage: x = ZZ['x'].gen()
            sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
            sage: P3 = HilbertPullback(HilbertModularGroup(K3))
            sage: P3.D() # abs tol 1e-10
            2.0121921726123237
        """
        n = self.group().base_ring().number_field().degree()
        return float(self.max_ideal_norm()**(1/n)*self._exp_matrix_BLambda_norm())

    #
    # A collection of bounds necessary for finding the closest cusp.
    #
    def _bound_for_closest_cusp(self):
        """
        This is the bound such that if a cusp is closer than this then it is the closest cusp.
        Reference: Lemma XX

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: P1 = HilbertPullback(HilbertModularGroup(5))
            sage: P1._bound_for_closest_cusp() # abs tol 1e-10
            0.19098300562505255
            sage: P2 = HilbertPullback(HilbertModularGroup(10))
            sage: P2._bound_for_closest_cusp() # abs tol 1e-10
            0.004389006498287334
            sage: x = ZZ['x'].gen()
            sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
            sage: P3 = HilbertPullback(HilbertModularGroup(K3))
            sage: P3._bound_for_closest_cusp() # abs tol 1e-10
            0.00532645515454193

        """
        n = self.group().base_ring().number_field().degree()
        return self.max_ideal_norm()**(-1) * 2**(-n/2.) * \
                                self._exp_matrix_BLambda_norm()**(-2*n)

    def _bound_for_sigma_embeddings(self, z, prec=16):
        """
        Bound for the embeddings of the denominator of the closest cusp to z
        Reference: Corollary 5.

        INPUT::

        - ``z`` -- point in the upper half-plane
        - ``prec`` -- the number of bits precision in the returned values.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_sigma_embeddings(z)
            [1.306, 1.377]

            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_sigma_embeddings(z)
            [4.415, 4.654]

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(1,1),CC(0,3)])
            sage: P3._bound_for_sigma_embeddings(z)
            [43.49, 43.49,  25.11]


        """
        self._check_upper_half_plane_element(z)
        n = self.group().base_ring().degree()
        d = self.D()*z.imag_norm()**(-1/(2*n))
        return [upper(d*y**(-1/2),prec=prec) for y in z.imag()]


    def _bound_for_sigma_coordinates(self, z, prec=16):
        """
        Bound `c` for the coordinates, with respect to the ring of integers, of the closest cusp to z
        Reference: Lemma XXX

        INPUT::

        - ``z`` -- point in the upper half-plane
        - ``prec`` -- the number of bits precision in the returned values.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_sigma_coordinates(z)
            [1.358, 1.200]

            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_sigma_coordinates(z)
             [4.534, 1.434]

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(1,1),CC(0,3)])
            sage: P3._bound_for_sigma_coordinates(z)
            [37.62, 7.227, 6.488]


        """
        self._check_upper_half_plane_element(z)
        n = self.group().base_ring().degree()
        d = self.D()*z.imag_norm()**(-1/(2*n))
        bounds = []
        B = self.basis_matrix_ideal().inverse()
        for i in range(n):
            bd = 0
            for j, y in enumerate(z.imag()):
                bd += B[i,j].abs()*y**(-1/2)
            bounds.append(upper(bd*d,prec=prec))
        return bounds

    def _bound_for_rho_embeddings(self,z, sigma, prec=16):
        """
        Bound for the embeddings of the numerator of the closest cusp to z.
        Reference: Corollary 5

            rho_i in [sigma_i x_i - delta*y_i**0.5 , sigma_i x_i - delta*y_i**0.5]
        INPUT::

        - ``z`` -- point in the upper half-plane
        - ``sigma`` -- list of complex embeddings of algebraic integer
        - ``prec`` -- the number of bits precision in the returned values.


        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_rho_embeddings(z, 0)
            [(-1.306, 1.306), (-1.239, 1.239)]
            sage: b1,b2=H1.base_ring().gens();
            sage: P1._bound_for_rho_embeddings(z, b1)
            [(-1.925, 0.6880), (0.3790, 2.857)]
            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_rho_embeddings(z, 0)
            [(-4.415, 4.415), (-4.188, 4.188)]
            sage: b1,b2=H2.base_ring().gens();
            sage: sage: P2._bound_for_rho_embeddings(z, b2)
            sage: P2._bound_for_rho_embeddings(z, b2)
            [(-7.577, 1.253), (-1.026, 7.351)]


            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(1,1),CC(0,3)])
            sage: P3._bound_for_rho_embeddings(z, 0)
             [(-43.49, 43.49), (-43.49, 43.49), (-75.32, 75.32)]
            sage: b1,b2,b3=H3.base_ring().gens()
            sage: P3._bound_for_rho_embeddings(z, b1)
             [(-43.49, 43.49), (-43.17, 43.81), (-75.32, 75.32)]

        """
        self._check_upper_half_plane_element(z)
        n = self.group().base_ring().degree()
        d = self.D()*z.imag_norm()**(-1.0/(2*n))
        if not isinstance(sigma,list):
            sigma = self.group().base_ring().number_field()(sigma)
            sigma = sigma.complex_embeddings()
        factor = 1.01
        bounds = []
        for i in range(n):
            y = z.imag()[i]
            xs = z.real()[i]*sigma[i]
            dy = d*y**(0.5)
            b0 = xs - dy
            b1 = xs + dy
            bounds.append((b0,b1))
        res = []
        for b0,b1 in bounds:
            # We bound the lower bound differently depending on whether it is positive or negative
            b0 = lower(b0,prec=prec)
            b1 = upper(b1,prec=prec)
            res.append((b0,b1))
        return res

    def _bound_for_rho_coordinates(self, z, sigma, prec=16):
        """
        Bound for the coordinates, with respect to the ring of integers, of the numerator of the closest cusp to z.
        Reference: Lemma XXX
            rho_i in [sigma_i x_i - delta*y_i**0.5 , sigma_i x_i - delta*y_i**0.5]

        INPUT::

        - ``z`` -- point in the upper half-plane
        - ``sigma`` -- complex embeddings of algebraic integer
        - ``prec`` -- the number of bits precision in the returned values.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P1._bound_for_rho_coordinates(z, 0)
            [2.571, 2.571]
            sage: b1,b2=H1.base_ring().gens();
            sage: P1._bound_for_rho_coordinates(z, b1)
            [5.839, 4.205]
            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.9)])
            sage: P2._bound_for_rho_coordinates(z, 0)
            [8.689, 8.689]
            sage: b1,b2=H2.base_ring().gens();
            sage: P2._bound_for_rho_coordinates(z, b2)
            [21.99, 21.99]


            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(1,1),CC(0,3)])
            sage: P3._bound_for_rho_coordinates(z, 0)
            [164.0, 164.0, 164.0]
            sage: b1,b2,b3=H3.base_ring().gens()
            sage: P3._bound_for_rho_coordinates(z, b1)
            [165.9, 164.0, 165.9]

        """
        self._check_upper_half_plane_element(z)
        n = self.group().base_ring().degree()
        d = self.D()*z.imag_norm()**(-1.0/(2*n))
        if not isinstance(sigma,list):
            sigma = self.group().base_ring().number_field()(sigma)
            sigma = sigma.complex_embeddings()
        bounds = []
        factor = 1.01
        sz = [(sigma[j] * z.real()[j]).abs() for j in range(n)]
        dy = [d * z.imag()[j] ** 0.5 for j in range(n)]
        for i in range(n):
            bd = 0
            for j, y in enumerate(z.imag()):
                bd += dy[j] + sz[j] * self.basis_matrix_ideal()[i, j].abs()

            bounds.append(upper(bd*factor,prec=prec))
        return bounds

    def _candidate_integers_sigma(self, z, domain='polytope', return_polyhedron=False, ideal_basis=None,
                                  lattice_basis=None, sorted=True):
        """
        Compute a list of candidates for the denominator of the closest cusp.

        INPUT::

        - ``z`` -- element of type UpperHalfPlaneProductelement_class
        - ``a`` -- ideal or algebraic integer (default = 1). If an integer is given then the ideal is the principal ideal.
        - ``domain`` -- string: 'polytope' (default), 'boundingbox',
        - ``return_polyhedron`` -- boolean
        - `ideal_basis` -- list or =None,
        - `lattice_basis` -- list of lists corresponding to a numerical basie of the lattice corresponding ot the ring of integers.
        - `sorted` -- boolean -- True to return a list sorted by norm first and then lexicographically with respect to embeddings.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1._candidate_integers_sigma(z)[0:3]
            [0, -1, 1]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.7)])
            sage: P1._candidate_integers_sigma(z)[0:5]
            [0, -1, 1, 1/2*a + 1/2, -1/2*a - 1/2]

            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P2._candidate_integers_sigma(z)[0:9]
            [0, -1, 1, -2, 2, -3, 3, a, -a]

            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,4),CC(0,5),CC(0,4)])
            sage: P3._candidate_integers_sigma(z)[0:5]
            [0, -1, 1, -2, 2]
        """
        if return_polyhedron:
            if domain == 'boundingbox':
                bounds = self._bound_for_sigma_coordinates(z)
                B = None
            elif domain == 'preimage':
                bounds = self._bound_for_sigma_embeddings(z)
                B = None
            else:
                bounds = self._bound_for_sigma_embeddings(z)
                B = self.basis_matrix_ideal().inverse()
            return self.polytope_from_bounds(bounds,B)
        # Else we call efficient methods to find candidates.
        if not lattice_basis or not ideal_basis:
            lattice_basis, ideal_basis = self._get_lattice_and_ideal_basis()
        coordinate_bounds = self._bound_for_sigma_coordinates(z)
        embedding_bounds = self._bound_for_sigma_embeddings(z)
        coordinate_bounds = [(-b, b) for b in coordinate_bounds]
        embedding_bounds = [(-b, b) for b in embedding_bounds]
        sigma_candidates_coordinates = lattice_elements_in_box(lattice_basis,
                                                               embedding_bounds,
                                                               coordinate_bounds)
        sigma_candidates = coordinates_to_ideal_elements(sigma_candidates_coordinates,
                                                             ideal_basis)
        if sorted:
            def absort(val):
                return (sum([abs(x)**2 for x in val.complex_embeddings()]),) + tuple(val.complex_embeddings())
            sigma_candidates.sort(key=absort)
        return sigma_candidates

    def _get_lattice_and_ideal_basis(self, ideal=None):
        """
        Compute an integral basis for an ideal as well as a basis matrix for the associated lattice in R^n in the form of a nested list.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: P1._get_lattice_and_ideal_basis()
            ([[1.00000000000000, -1.61803398874989],
              [1.00000000000000, 0.618033988749895]],
             [1, 1/2*a - 1/2])
            sage: P1._get_lattice_and_ideal_basis(H1.base_ring().ideal(2))
            ([[2.00000000000000, -3.23606797749979], [2.00000000000000, 1.23606797749979]],
            [2, a - 1])
            sage: H2=HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: P2._get_lattice_and_ideal_basis()
             ([[1.00000000000000, -3.16227766016838], [1.00000000000000, 3.16227766016838]],
              [1, a])
            sage: a=H2.base_ring().ideal(H2.base_ring().gen(1)+1)
            sage: P2._get_lattice_and_ideal_basis(a)
             ([[9.00000000000000, -2.16227766016838], [9.00000000000000, 4.16227766016838]],
              [9, a + 1])
            sage: x = ZZ['x'].gen()
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: P3=HilbertPullback(H3)
            sage: P3._get_lattice_and_ideal_basis()
            ([[1.00000000000000, -5.98606258583498, 2.28229423189948],
              [1.00000000000000, -0.0277783731902446, -7.67566891172438],
              [1.00000000000000, 6.01384095902523, 6.39337467982490]],
             [1, a, 1/3*a^2 + 1/3*a - 23/3])
            sage: c = P3.group().base_ring().class_group()[1].ideal()
            sage: P3._get_lattice_and_ideal_basis(c)
            ([[2.00000000000000, -4.98606258583498, 3.28229423189948],
              [2.00000000000000, 0.972221626809755, -6.67566891172438],
              [2.00000000000000, 7.01384095902523, 7.39337467982490]],
             [2, a + 1, 1/3*a^2 + 1/3*a - 20/3])



        """
        if not ideal:
            ideal = self._construct_ideal(1)
        ideal_basis = ideal.basis()
        lattice_basis = self.basis_matrix_ideal(ideal)
        n = len(lattice_basis[0])
        # Make lattice basis to a nested list to avoid creation of FreeModule elements
        lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
        return lattice_basis, ideal_basis

    def _candidate_integers_rho(self, z, sigma, a=1, domain='polytope',return_polyhedron=False,
                                ideal_basis=None,
                                lattice_basis=None, sorted=True):
        """
        Compute a list of candidates for the denominator of the closest cusp.

         INPUT::

        - ``z`` -- element of type UpperHalfPlaneProductelement_class
        - ``sigma`` -- algebraic integer
        - ``a`` -- ideal or algebraic integer (default = 1). If an integer is given then the ideal is the principal ideal.
        - ``algorithm`` -- string (either 'coordinates' or 'embeddings')

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P1._candidate_integers_rho(z,1)
            [0, 1, 2]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,0.25)])
            sage: P1._candidate_integers_rho(z,1)
            [1, 1/2*a + 1/2, -1/2*a + 3/2]
            sage: H2 = HilbertModularGroup(10)
            sage: P2 = HilbertPullback(H2)
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(1,1)])
            sage: P2._candidate_integers_rho(z,1)
            [0, -1, 1, -2, 2, -3, 3, a, -a, a + 1, -a + 1, a + 2, -a + 2, 4, 5]

            # sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            # sage: P3=HilbertPullback(H3)
            # sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1),CC(0,1)])
            # sage: sigma = H3.base_ring().gens()[1]
            # sage: P3._candidate_integers_rho(z, sigma)


        """
        if return_polyhedron:
            if domain == 'boundingbox':
                bounds = self._bound_for_rho_coordinates(z, sigma)
                B = None
            elif domain == 'preimage':
                bounds = self._bound_for_rho_embeddings(z, sigma)
                B = None
            else:
                bounds = self._bound_for_rho_embeddings(z, sigma)
                B = self.basis_matrix_ideal().inverse()
            return self.polytope_from_bounds(bounds,B)
        if not lattice_basis or not ideal_basis:
            lattice_basis, ideal_basis = self._get_lattice_and_ideal_basis()
        # if not isinstance(sigma,list):
        #     sigma = self.group().base_ring().number_field()(sigma)
        #     sigma = sigma.complex_embeddings()
        rho_coordinate_bounds = self._bound_for_rho_coordinates(z, sigma)
        rho_coordinate_bounds = [(-b,b) for b in rho_coordinate_bounds]
        rho_embedding_bounds = self._bound_for_rho_embeddings(z, sigma)

        rho_candidates_coordinates = lattice_elements_in_box(lattice_basis,
                                                 rho_embedding_bounds,
                                                 rho_coordinate_bounds)
        rho_candidates = coordinates_to_ideal_elements(rho_candidates_coordinates,
                                                         ideal_basis)
        if sorted:
            def absort(val):
                return (sum([abs(x)**2 for x in val.complex_embeddings()]),) + tuple(val.complex_embeddings())
            rho_candidates.sort(key=absort)
        return rho_candidates

    def _candidate_closest_cusps(self, z, a=None, domain='polytope', as_cusps=False):
        r"""
        Find candidates for the closest cusp.

        INPUT::

        - ``z`` -- point in the upper half-plane
        - ``a`` -- ideal
        - `domain` -- string, either 'polytope' or 'boundingbox'
        - `as_cusps` -- boolean - deafult False, set to True to return a list of cusps instead of tuples.


        EXAMPLES::
            sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
            sage: H1 = HilbertModularGroup(5)
            sage: P1 = HilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
            sage: P1._candidate_closest_cusps(z)
            [(1, 0), (0, 1), (-1, -1), (1, -1)]

        """
        lattice_basis, ideal_basis = self._get_lattice_and_ideal_basis()
        cusp_candidates = [(1,0)]
        sigmas = self._candidate_integers_sigma(z,domain=domain,lattice_basis=lattice_basis,ideal_basis=ideal_basis)
        for sigma in sigmas:
            if sigma == 0:
                continue
            rhos = self._candidate_integers_rho(z, sigma, domain=domain,lattice_basis=lattice_basis,ideal_basis=ideal_basis)
            for rho in rhos:
                if rho == 0 and (0, 1) not in cusp_candidates:
                    cusp_candidates.append((0, 1))
                elif rho == 0:
                    continue
                if (-rho,-sigma) not in cusp_candidates:
                    cusp_candidates.append((rho, sigma))
        if as_cusps:
            # filter away doubles
            cusps = []
            K = self.group().base_ring().number_field()
            for rho, sigma in cusp_candidates:
                c = NFCusp(K,rho,sigma)
                if c not in cusps:
                    cusps.append(c)
            return cusps
        else:
            return cusp_candidates

# H1=HilbertModularGroup(5)
# P1=HilbertPullback(H1)
# H2=HilbertModularGroup(10)
# P2=HilbertPullback(H2)
# from sage.all import NumberField, QQ, RDF
#
# x = QQ['x'].gen()
# H3=HilbertModularGroup(NumberField(x**3-36*x-1, names='a'))
# P3=HilbertPullback(H3)
# from sage.all import ceil
# def integral_coordinates_in_box(bounds):
#     """
#     Find all integral points inside a box in R^n.
#
#     INPUT:
#     - `bounds` -- a list of upper and lower bounds defining a parallelepiped of positive volume.
#
#     """
#     if not isinstance(bounds, list):
#         raise ValueError
#     n = len(bounds)
#     coordinates = []
#     for i in range(n):
#         lower = ceil(bounds[i][0])
#         upper = floor(bounds[i][1])
#         if lower >= upper:
#             raise ValueError("Bounds must give interval of positive length.")
#         coordinates.append(range(lower, upper + 1))
#     # result = [tuple(x) for x in cartesian_product(coordinates)]
#     # return result
#     return list(cartesian_product(coordinates))
#
# def lattice_elements_in_box(lattice_basis, lattice_bounds, coordinate_bounds):
#         # bounds_coordinates, bounds_embeddings,
#         #                   basis, basis_embeddings):
#     """
#     Find all coordinates of elements in a lattice in R^n restricted to a specific box and with
#     coordinates in another box.
#
#     """
#     coordinates = integral_coordinates_in_box(coordinate_bounds)
#     result = []
#     n = len(lattice_basis[0])
#     for coordinate_vector in coordinates:
#         is_within_bounds = True
#         for i in range(n):
#             alpha_i = 0.0
#             for j in range(n):
#                 alpha_i = alpha_i + lattice_basis[i][j] * coordinate_vector[j]
#             if alpha_i < lattice_bounds[i][0] or alpha_i > lattice_bounds[i][1]:
#                 # We need to discard this
#                 is_within_bounds = False
#                 break
#         # If we are within the bounds we add the number field element.
#         if is_within_bounds:
#             result.append(coordinate_vector)
#     return result
#
# def coordinates_to_ideal_elements(coordinates,ideal_basis):
#     result = []
#     for coordinate_vector in coordinates:
#         element = 0
#         for i,b in enumerate(ideal_basis):
#             element += b * coordinate_vector[i]
#         result.append(element)
#     return result
#     # return [
#     #     sum(b * coordinate_vector[i] for i,b in enumerate(ideal_basis))
#     #     for coordinate_vector in coordinates]
#
# def find_sigma_candidates(p, z, ideala):
#     ideal_basis = ideala.integral_basis()
#     lattice_basis = p.basis_matrix_ideal()
#     n = len(lattice_basis[0])
#     # Make lattice basis to a nested list to avoid creation of FreeModule elements
#     lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
#     coordinate_bounds = p._bound_for_sigma_coordinates(z)
#     embedding_bounds = p._bound_for_sigma_embeddings(z)
#     coordinate_bounds = [(-b,b) for b in coordinate_bounds]
#     embedding_bounds = [(-b,b) for b in embedding_bounds]
#     sigma_candidates_coordinates = lattice_elements_in_box(lattice_basis,
#                                                embedding_bounds,
#                                                coordinate_bounds)
#     sigma_candidates = coordinates_to_ideal_elements(sigma_candidates_coordinates,
#                                                                  ideal_basis)
#     return sigma_candidates
#
# def find_candidate_cusps(p, z,ideala):
#     ideal_basis = ideala.integral_basis()
#     lattice_basis = p.basis_matrix_ideal()
#     n = len(lattice_basis[0])
#     # Make lattice basis to a nested list to avoid creation of FreeModule elements
#     lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
#     coordinate_bounds = p._bound_for_sigma_coordinates(z)
#     embedding_bounds = p._bound_for_sigma_embeddings(z)
#     coordinate_bounds = [(-b,b) for b in coordinate_bounds]
#     embedding_bounds = [(-b,b) for b in embedding_bounds]
#     sigma_candidates_coordinates = lattice_elements_in_box(lattice_basis,
#                                                embedding_bounds,
#                                                coordinate_bounds)
#     sigma_candidates = coordinates_to_ideal_elements(sigma_candidates_coordinates,
#                                                              ideal_basis)
#     result = [(1,0)]
#     # print("sigma candidates=",len(sigma_candidates))
#     for s in sigma_candidates:
#         if s == 0:
#             continue
#         rho_coordinate_bounds = p._bound_for_rho_coordinates(z, s)
#         rho_coordinate_bounds = [(-b,b) for b in rho_coordinate_bounds]
#         rho_embedding_bounds = p._bound_for_rho_embeddings(z, s)
#         # rho_embedding_bounds = [(-b, b) for b in rho_embedding_bounds]
#
#         rho_candidates_coordinates = lattice_elements_in_box(lattice_basis,
#                                                  rho_embedding_bounds,
#                                                  rho_coordinate_bounds)
#         rho_candidates = coordinates_to_ideal_elements(rho_candidates_coordinates,
#                                                          ideal_basis)
#
#         for r in rho_candidates:
#             #c = NFCusp(p.number_field(),r,s)
#             if (r,s) in result:
#                 continue
#             if (-r,-s) in result:
#                 continue
#             # if c in result:
#             #     continue
#             if r == 0 and ((0,1) in result or s != 1):
#                continue
#             result.append((r,s))
#     return result
#
# def find_closest_cusp(p, z, ideala, return_multiple=False):
#     cusp_candidates = find_candidate_cusps(p, z, ideala)
#     min_cusp = cusp_candidates[0]
#     min_d = distance_to_cusp(p,min_cusp[0],min_cusp[1], z)
#     if return_multiple:
#         min_cusp = [min_cusp]
#     for c in cusp_candidates[1:]:
#         d = distance_to_cusp(p,c[0],c[1], z)
#         if d < min_d:
#             if return_multiple:
#                 min_cusp = [c]
#             else:
#                 min_cusp = c
#         if d == min_d and return_multiple:
#             min_cusp.append(c)
#     return min_cusp
#
# def distance_to_cusp(p,r,s,z):
#     ideal_rs = p.number_field().ideal(r,s)
#     x = z.real()
#     y = z.imag()
#     n = len(y)
#     if hasattr(r,'complex_embeddings'):
#         r = r.complex_embeddings()
#     else:
#         r= [r]*n
#     if hasattr(s,'complex_embeddings'):
#         s = s.complex_embeddings()
#     else:
#         s = [s]*n
#     d = 1
#     for i in range(n):
#         d = d* ((x[i] * s[i] + r[i]) ** 2 * y[i] ** -1 + s[i] ** 2 * y[i])
#     d = ideal_rs.norm()**-1*d.sqrt()
#     return d