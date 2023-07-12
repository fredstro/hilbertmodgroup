r"""
Classes for (projective) Hilbert modular groups PSL_2(OK)

EXAMPLES::

    sage: from hilbert_modgroup.all import HilbertModularGroup
    sage: HilbertModularGroup(5)
    Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?
    sage: HilbertModularGroup(QuadraticField(5))
    Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?


AUTHORS:
- Fredrik Stromberg (2021)


"""
import sage
from sage.categories.groups import Groups
from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic
from sage.modular.cusps_nf import NFCusp
from sage.rings.infinity import infinity
from sage.rings.number_field.number_field import QuadraticField, CyclotomicField
from sage.all import latex, Integer, Matrix, matrix
from sage.misc.cachefunc import cached_method
from sage.rings.number_field.order import is_NumberFieldOrder

# from .upper_half_plane import ComplexPlaneOtimesK
from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement__class
from .hilbert_modular_group_element import HilbertModularGroupElement

import logging

logger = logging.getLogger(__name__)
logger.setLevel(10)


def is_HilbertModularGroup(x):
    """
    Return `True` if ``x`` is an instance of a HilbertModularGroup

    INPUT:

    - ``x`` -- something to test if it is a Hilbert modular group or not

    OUTPUT:

    - boolean

    EXAMPLES::

        sage: from hilbert_modgroup.all import HilbertModularGroup,is_HilbertModularGroup
        sage: is_HilbertModularGroup(1)
        False
        sage: H = HilbertModularGroup(5)
        sage: is_HilbertModularGroup(H)
        True
    """
    return isinstance(x, HilbertModularGroup_class)


def HilbertModularGroup(number_field, projective=True):
    r"""
        Create the Hilbert modular group over the ring of integers in the given number field


        INPUT:

        - ``number_field`` (NumberField) -- a totally real number field or positive integer.
                                            If a positive integer D is specified
                                            then the number field $Q(\sqrt(D))$ is used.
        - ``projective`` (bool) - True if you want PSL(2,K) and False for SL(2,K)


        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5)
            Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?
            sage: HilbertModularGroup(QuadraticField(5))
            Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?


        TESTS::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5)
            Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?
            sage: HilbertModularGroup(QuadraticField(5))
            Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?


    """
    if isinstance(number_field, (int, Integer)) and number_field > 0:
        ring = QuadraticField(number_field).ring_of_integers()
    elif isinstance(number_field, sage.rings.number_field.number_field_base.NumberField) \
            and number_field.is_totally_real():
        ring = number_field.ring_of_integers()
    else:
        raise ValueError("The input must be a totally real Number Field or a positive integer")
    if not projective:
        raise NotImplementedError("Only PSL2 is implemented at the moment.")
    degree = Integer(2)
    name = f'Hilbert Modular Group PSL({degree}) over {ring}'
    ltx = 'PSL({0}, {1})'.format(degree, latex(ring))
    return HilbertModularGroup_class(base_ring=ring, sage_name=name, latex_string=ltx)


class HilbertModularGroup_class(LinearMatrixGroup_generic):
    r"""
    Class for Hilbert modular groups, here defined as either PSL(2) (default) or  SL(2)
    over rings of integers in totally real number fields.


    """

    Element = HilbertModularGroupElement

    def __init__(self, base_ring, sage_name, latex_string):
        r"""
         Init a Hilbert modular group over the ring of integers in the given number field


        INPUT:
        - ``base_ring`` - ring
        - ``sage_name`` - string
        - ``latex_string`` - string

        EXAMPLES::
            sage: from hilbert_modgroup.hilbert_modular_group_class import *
            sage: OK=QuadraticField(5).OK()
            sage: name = f'Hilbert Modular Group PSL(2) over {OK}'
            sage: ltx = f'PSL(2, {latex(OK)})'
            sage: HilbertModularGroup_class(base_ring=OK,sage_name=name,latex_string=ltx)
            Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?
            sage: H1=HilbertModularGroup(5)
            sage: TestSuite(H1).run()
            sage: H1(1)
            [1 0]
            [0 1]
            sage: H1(2)
            Traceback (most recent call last):
            ...
            TypeError: matrix must have determinant 1
            sage: H1([1,1,0,1])
            [1 1]
            [0 1]
            sage: H1([1,H1.base_ring().gens()[0],0,1])
            [          1 1/2*a + 1/2]
            [          0           1]
        """
        if not is_NumberFieldOrder(base_ring) or not base_ring.number_field().is_totally_real():
            raise ValueError("Input (={0}) can not be used to create a Hilbert modular group. " +
                             "Need an order of a totally real number field")
        # Instance data related to elliptic elements
        self._elliptic_elements_traces = []
        self._elliptic_elements_orders = []
        self._elliptic_elements_traces_of_orders = {}
        self._elliptic_elements_fields_of_orders = {}
        self._elliptic_elements_orders_of_traces = {}
        # Instance data related to cusps
        self._ncusps = None
        self._cusps = []
        self._ideal_cusp_representatives = []
        self._cusp_normalizing_maps = {}
        self._cusp_normalizing_maps_inverse = {}
        # At the moment we only deal with full level (1)
        self._level = base_ring.fractional_ideal(1)
        super(HilbertModularGroup_class, self).__init__(degree=Integer(2), base_ring=base_ring,
                                                        special=True,
                                                        sage_name=sage_name,
                                                        latex_string=latex_string,
                                                        category=Groups().Infinite(),
                                                        invariant_form=None)

    @cached_method
    def generators(self, algorithm='standard'):
        r"""
        Return a list of generators of self.

        INPUT:


        - ``algorithm`` (string) either 'standard' or 'elementary'.
            If 'elementary' is given return a set of generators
         consisting of elementary (i.e. upper- and lower-triangular) matrices.
         Otherwise return a set of reflections and translations.


        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: H.generators()
            [
            [ 0 -1]  [          1 1/2*a + 1/2]  [1 a]
            [ 1  0], [          0           1], [0 1]
            ]

            sage: H.generators(algorithm='elementary')
            [
            [          1 1/2*a + 1/2]  [          1           0]  [1 a]  [1 0]
            [          0           1], [1/2*a + 1/2           1], [0 1], [a 1]
            ]


        """
        if algorithm == 'standard':
            gens = [self.S()]
            for x in self.base_ring().basis():
                gens.append(self.T(x))
        elif algorithm == 'elementary':
            gens = []
            for x in self.base_ring().basis():
                gens.append(self.T(x))
                gens.append(self.L(x))
        else:
            raise ValueError("Unknown algorithm '{0}'. Expected one of 'standard' or 'elementary'")
        return gens

    @cached_method
    def S(self):
        """
        Return the element S = ( 0 & -1 // 1 & 0 ) of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5).S()
            [ 0 -1]
            [ 1  0]
            sage: HilbertModularGroup(10).S()
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

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: H.T()
            [1 1]
            [0 1]
            sage: u0,u1=H.base_ring().number_field().unit_group().gens()
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

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H = HilbertModularGroup(5)
            sage: H.L(1)
            [1 0]
            [1 1]
            sage: u0,u1=H.base_ring().number_field().unit_group().gens()
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
        Return the element U=( u & 0 // 0 & u**-1 ) of self.

        INPUT:
        - `u` unit in self.base_ring()

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H = HilbertModularGroup(5)
            sage: H.E(1)
            [1 0]
            [0 1]
            sage: u0,u1=H.base_ring().number_field().unit_group().gens()
            sage: H.E(u0)
            [-1  0]
            [ 0 -1]
            sage: H.E(u0*u1)
            [1/2*a - 1/2           0]
            [          0 1/2*a + 1/2]

        """
        return self([u, 0, 0, u ** -1])

    def gens(self, algorithm='standard'):
        r"""
            Return a tuple of generators for this Hilbert modular group.

            The generators need not be minimal. For arguments, see :meth:`~generators`.

        INPUT:

        - ``algorithm`` -- string (default='standard') give the algorithm to compute the generators

        NOTE: Different 'algorithms' give different choices of generators.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5).gens()
            (
            [ 0 -1]  [          1 1/2*a + 1/2]  [1 a]
            [ 1  0], [          0           1], [0 1]
            )

        """
        return tuple(self.generators(algorithm))

    def ngens(self, algorithm='standard'):
        r"""
            Return the number of generators of self as given by the function 'gens'.

        INPUT:

        - ``algorithm`` -- string (default='standard') give the algorithm to compute the generators

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5).ngens()
            3

        """
        return len(self.generators(algorithm))

    def gen(self, i):
        r"""
        Return the i-th generator of self, i.e. the i-th element of the
        tuple self.gens().

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5).gen(1)
            [          1 1/2*a + 1/2]
            [          0           1]
        """
        return self.generators()[i]

    def random_element(self, *args, **kwds):
        r"""
        Return a 'random' element of this Hilbert Modular Group.

        INPUT:

        - `args`, `kwds` -- arguments passed to the base ring's random element function
                            and are in turn passed to the random integer function.
                            See the documentation for "ZZ.random_element()" for details.


        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H = HilbertModularGroup(5)
            sage: A = H.random_element()
            sage: A in H
            True

        """
        a = self.base_ring().random_element(**kwds)
        b = self.base_ring().random_element(**kwds)
        return self.T(a)*self.L(b)

    def level(self):
        """
        Return the level of this Hilbert modular group (currently only (1))

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: HilbertModularGroup(5).level()
            Fractional ideal (1)

        """
        return self._level

    def ncusps(self):
        """
        Return number of cusps of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.ncusps()
            1

            sage: H2=HilbertModularGroup(10)
            sage: H2.ncusps()
            2

            sage: var('x')
            x
            sage: K = NumberField(x^3-36*x-1, names='a')
            sage: H3=HilbertModularGroup(K)
            sage: H3.ncusps()
            5

            sage: K4 = NumberField(x**4 - 17*x**2 + 36,names='a'); a=K4.gen()
            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.ncusps()
            2
        """
        if not self._ncusps:
            self._ncusps = self.base_ring().class_number()
        return self._ncusps

    @cached_method
    def cusps(self):
        """
        A set of cusp representatives of self.


        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.cusps()
            [Cusp Infinity of Number Field ... polynomial x^2 - 5 with a = 2.236067977499790?]

            sage: H2=HilbertModularGroup(10)
            sage: H2.cusps()
            [Cusp Infinity of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?,
            Cusp [2: a] of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?]
            sage: x=var('x')
            sage: K = NumberField(x^3-36*x-1, names='a')
            sage: H3=HilbertModularGroup(K); H3
            Hilbert Modular Group ... x^3 - 36*x - 1
            sage: H3.cusps()
             [Cusp Infinity of Number Field in a with defining polynomial x^3 - 36*x - 1,
             Cusp [2: a + 1] of Number Field in a with defining polynomial x^3 - 36*x - 1,
             Cusp [3: 1/3*a^2 + 1/3*a - 26/3] of Number Field ... polynomial x^3 - 36*x - 1,
             Cusp [2: 1/3*a^2 + 1/3*a - 23/3] of Number Field ... polynomial x^3 - 36*x - 1,
             Cusp [6: 1/3*a^2 + 1/3*a - 26/3] of Number Field ... polynomial x^3 - 36*x - 1]

            sage: K4 = NumberField(x**4 - 17*x**2 + 36,names='a'); a=K4.gen()
            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.cusps()
            [Cusp Infinity of Number Field in a with defining polynomial x^4 - 17*x^2 + 36,
             Cusp [2: a + 1] of Number Field in a with defining polynomial x^4 - 17*x^2 + 36]


        """
        for a in self.ideal_cusp_representatives():
            logger.debug("Set cusp info for ideal a={0}".format(a))
            if a.is_trivial():
                ca = NFCusp(self.base_ring().number_field(),
                            self.base_ring()(1),
                            self.base_ring()(0),
                            lreps=self.ideal_cusp_representatives())
            else:
                ag = a.gens_reduced()
                ca = NFCusp(self.base_ring().number_field(), ag[0], ag[1],
                            lreps=self.ideal_cusp_representatives())
            self._cusps.append(ca)
            if ca.ideal() != a:
                raise ArithmeticError("Failed to associate a cusp to ideal {0}".format(a))
        return self._cusps

    def ideal_cusp_representatives(self):
        r"""
        Return a list of ideals corresponding to cusp representatives, i.e.
        ideal representatives of ideal classes.

        Note: We choose an ideal of smallest norm in each class.
            If the ideal given by sage is already minimal we return this.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.ideal_cusp_representatives()
            [Fractional ideal (1)]

            sage: H2=HilbertModularGroup(10)
            sage: H2.ideal_cusp_representatives()
            [Fractional ideal (1), Fractional ideal (2, a)]
            sage: x=var('x')
            sage: K = NumberField(x^3-36*x-1, names='a')
            sage: H3=HilbertModularGroup(K); H3
            Hilbert Modular Group ... x^3 - 36*x - 1
            sage: H3.ideal_cusp_representatives()
            [Fractional ideal (1),
             Fractional ideal (2, a + 1),
             Fractional ideal (3, 1/3*a^2 + 1/3*a - 26/3),
             Fractional ideal (2, 1/3*a^2 + 1/3*a - 23/3),
             Fractional ideal (6, 1/3*a^2 + 1/3*a - 26/3)]
            sage: K4 = NumberField(x**4 - 17*x**2 + 36,names='a'); a=K4.gen()
            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.ideal_cusp_representatives()
            [Fractional ideal (1), Fractional ideal (2, a + 1)]


        """
        if not self._ideal_cusp_representatives:
            self._ideal_cusp_representatives = []

            def _find_equivalent_ideal_of_minimal_norm(c):
                for a in self.base_ring().number_field().ideals_of_bdd_norm(c.norm() - 1).items():
                    for ideala in a[1]:
                        if (ideala * c ** -1).is_principal():
                            if c.norm() <= ideala.norm():
                                return c
                            return ideala
                return c

            for ideal_class in self.base_ring().class_group():
                c = ideal_class.ideal().reduce_equiv()
                # NOTE: Even though we use 'reduce_equiv' we are not guaranteed a representative
                #       with minimal **norm**
                #       To make certain we choose a representative of minimal norm explicitly.
                c = _find_equivalent_ideal_of_minimal_norm(c)
                self._ideal_cusp_representatives.append(c)
            # We finally sort all representatives according to norm.
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

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: c = NFCusp(H1.base_ring().number_field(),2,4)
            sage: H1.cusp_representative(c)
            Cusp Infinity of Number Field ... polynomial x^2 - 5 with a = 2.236067977499790?

            sage: H2=HilbertModularGroup(10)
            sage: c = NFCusp(H2.base_ring().number_field(),2,4)
            sage: H2.cusp_representative(c)
            Cusp Infinity of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?

            sage: a = H2.base_ring().number_field().gen()
            sage: x,y = 3*a - 10, a - 4
            sage: c = NFCusp(H2.base_ring().number_field(),x,y)
            sage: H2.cusp_representative(c)
            Cusp [2: a] of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?

            sage: x = ZZ['x'].gen()
            sage: K = NumberField(x^3-36*x-1, names='a')
            sage: H3=HilbertModularGroup(K)
            sage: a = K.gen()
            sage: c = NFCusp(H3.base_ring().number_field(),2,3)
            sage: H3.cusp_representative(c)
            Cusp Infinity of Number Field in a with defining polynomial x^3 - 36*x - 1
            sage: x,y = 3*a - 10, a - 4
            sage: c = NFCusp(H3.base_ring().number_field(),16,a+3)
            sage: H3.cusp_representative(c)
            Cusp [2: 1/3*a^2 + 1/3*a - 23/3] of Number Field ... polynomial x^3 - 36*x - 1


        """
        for c in self.cusps():
            if return_map:
                t, B = cusp.is_Gamma0_equivalent(c, self.level(), Transformation=True)
                if t:
                    return c, self(B)
            elif cusp.is_Gamma0_equivalent(c, self.level()):
                return c
        raise ArithmeticError(f"Could not find cusp representative for {cusp}")

    # Functions for elliptic elements

    def _compute_traces_of_elliptic_elements(self):
        r"""
        Compute all possible traces of elliptic elements for self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1._compute_traces_of_elliptic_elements()
            [-1, 0, -1/2*a - 1/2, 1/2*a - 1/2, 1, -1/2*a + 1/2, 1/2*a + 1/2]
            sage: H2=HilbertModularGroup(10)
            sage: H2._compute_traces_of_elliptic_elements()
            [-1, 0, 1]


        ALGORITHM: (from the paper "Dimension formulas...", by Boylan, Skoruppa and Stromberg)

        If m is an elliptic element then:
        1. The order, l, must satisfy euler_phi(l)=2d where  d | n = degree of the number field.
        2. The trace t must satisfy z+z^-1 where z is an l-th root of unity.
        3. The subfield QQ(z+z^-1) of QQ(z) must be a subfield of the base field K
        Conversely, if t is such then there is an elliptic element with t as trace.
        These two conditions therefore characterise the traces of elliptic points completely.

        """
        from sage.arith.misc import euler_phi
        if self._elliptic_elements_traces:
            return self._elliptic_elements_traces
        K = self.base_ring().number_field()
        n = K.degree()
        list_of_possible_orders = [o for o in range(2, 8*n**2+1) if euler_phi(o).divides(2*n)]
        for o in list_of_possible_orders:
            C = CyclotomicField(o, 'z')
            z = C.gen()
            F, emb = C.subfield(z+z**-1, 't')
            if not F.embeddings(K):
                continue
            t = F.gen()
            possible_traces_in_K = [s(t) for s in F.embeddings(K)]
            logger.debug("F={0}".format(F))
            logger.debug("z={0}".format(z))
            logger.debug("t={0}".format(t))
            logger.debug("|F.emb(K)|={0}".format(len(F.embeddings(K))))
            logger.debug("t={0}".format(possible_traces_in_K))
            traces_of_order_o = []
            for st in possible_traces_in_K:
                # Make sure that the trace have all embeddings with absolute value <2
                test = [x for x in st.complex_embeddings() if abs(x) >= 2]
                logger.debug("st={0}".format(st))
                logger.debug("test={0}".format(test))
                if test:
                    continue
                # We want to choose a representative since the point only depends
                # on t^2 we don't need t and -t
                if st not in traces_of_order_o:
                    traces_of_order_o.append(st)
            if traces_of_order_o:
                self._elliptic_elements_traces.extend(traces_of_order_o)
                self._elliptic_elements_fields_of_orders[o] = F
                self._elliptic_elements_traces_of_orders[o] = traces_of_order_o
                self._elliptic_elements_orders.append(o)
        self._elliptic_elements_orders_of_traces = {
            value: key
            for key in self._elliptic_elements_traces_of_orders
            for value in self._elliptic_elements_traces_of_orders[key]}
        return self._elliptic_elements_traces

    def orders_of_elliptic_elements(self):
        r"""
        Return a list of possible orders of elliptic elements.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.orders_of_elliptic_elements()
            [3, 4, 5, 6, 10]

            sage: H2=HilbertModularGroup(10)
            sage: H2.orders_of_elliptic_elements()
            [3, 4, 6]

            sage: var('x')
            x
            sage: H4=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: H4.orders_of_elliptic_elements()
            [3, 4, 6]

            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.orders_of_elliptic_elements()
            [3, 4, 5, 6, 10]

        """
        if not self._elliptic_elements_orders:
            self._compute_traces_of_elliptic_elements()
        return self._elliptic_elements_orders

    def traces_of_elliptic_elements(self):
        r"""
        Return a list of possible traces of elliptic elements.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: H.traces_of_elliptic_elements()
            [-1, 0, -1/2*a - 1/2, 1/2*a - 1/2, 1, -1/2*a + 1/2, 1/2*a + 1/2]

            sage: H2=HilbertModularGroup(10)
            sage: H2.traces_of_elliptic_elements()
            [-1, 0, 1]

            sage: var('x')
            x
            sage: K = NumberField(x^3-36*x-1, names='a')
            sage: H=HilbertModularGroup(K)
            sage: H2.traces_of_elliptic_elements()
            [-1, 0, 1]

            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.orders_of_elliptic_elements()
            [3, 4, 5, 6, 10]
            sage: H4.traces_of_elliptic_elements()
            [-1,
             0,
             1/12*a^3 - 11/12*a - 1/2,
             -1/12*a^3 + 11/12*a - 1/2,
             1,
             1/12*a^3 - 11/12*a + 1/2,
             -1/12*a^3 + 11/12*a + 1/2]

        """
        if not self._elliptic_elements_traces:
            self._compute_traces_of_elliptic_elements()
        return self._elliptic_elements_traces

    def traces_of_elliptic_elements_of_order(self, o):
        r"""
        Return a list of traces of elliptic elements of given order.

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.traces_of_elliptic_elements_of_order(4)
            [0]
            sage: H1.traces_of_elliptic_elements_of_order(7)
            []

            sage: H2=HilbertModularGroup(10)
            sage: H2.traces_of_elliptic_elements_of_order(3)
            [-1]

            sage: var('x')
            x
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: H3.traces_of_elliptic_elements_of_order(6)
            [1]

            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.traces_of_elliptic_elements_of_order(10)
            [1/12*a^3 - 11/12*a + 1/2, -1/12*a^3 + 11/12*a + 1/2]
            sage: H4.traces_of_elliptic_elements_of_order(5)
            [1/12*a^3 - 11/12*a - 1/2, -1/12*a^3 + 11/12*a - 1/2]

        """
        if not self._elliptic_elements_traces_of_orders:
            self._compute_traces_of_elliptic_elements()
        return self._elliptic_elements_traces_of_orders.get(o, [])

    def order_of_elliptic_element_of_trace(self, t):
        r"""
        Return the order of elliptic elements of a given trace.
        Returns None if no elliptic element with this trace exists.

        INPUT:
        - `t` number field element

        OUTPUT:
        - integer or None



        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.order_of_elliptic_element_of_trace(0)
            4
            sage: H1.order_of_elliptic_element_of_trace(2)


            sage: H2=HilbertModularGroup(10)
            sage: H2.order_of_elliptic_element_of_trace(-1)
            3

            sage: var('x')
            x
            sage: H3=HilbertModularGroup(NumberField(x^3-36*x-1, names='a'))
            sage: H2.order_of_elliptic_element_of_trace(1)
            6

            sage: K4 = NumberField(x**4 - 17*x**2 + 36,names='a'); a=K4.gen()
            sage: H4=HilbertModularGroup(NumberField(x**4 - 17*x**2 + 36,names='a'))
            sage: H4.order_of_elliptic_element_of_trace(1)
            6
            sage: H4.order_of_elliptic_element_of_trace(1/12*a^3 - 11/12*a + 1/2)
            10

        """
        if not self._elliptic_elements_orders_of_traces:
            self._compute_traces_of_elliptic_elements()
        return self._elliptic_elements_orders_of_traces.get(t, None)

    #
    # Functions for working with cusps.
    #

    def cusp_normalizing_map(self, cusp, inverse=False, check=False):
        r"""
        Given a cusp (a:c) Return a matrix A = [[ a ,b ], [c , d]] in SL(2,K) such that
        A(Infinity)=(a:c) and b, d in self.base_ring().ideal(a,c)**-1

        INPUT:

        - ``cusp`` -- Instance of NFCusp
        - ``inverse`` -- bool (default: False) set to True to return the inverse map
        - ``check`` -- bool (default: False) set to True to check the result

        NOTE: The sage function NFCusp.ABmatrix() returns a matrix with determinant which is not
            necessarily equal to 1 even though 1 is a generator of the ideal (1)=(a,c)*(a,c)**-1

        If inverse = True then return A^-1

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H1=HilbertModularGroup(5)
            sage: H1.cusp_normalizing_map(H1.cusps()[0])
            [1 0]
            [0 1]
            sage: H1.cusp_normalizing_map(NFCusp(H1.base_ring().number_field(),0,1))
            [ 0 -1]
            [ 1  0]
            sage: H1.cusp_normalizing_map(NFCusp(H1.base_ring().number_field(),1,1))
            [ 1 -1]
            [ 1  0]
            sage: a=H1.base_ring().number_field().gen()
            sage: H1.cusp_normalizing_map(NFCusp(H1.base_ring().number_field(),a,1+a))
            [    a a - 1]
            [a + 1     a]

            sage: H2=HilbertModularGroup(10)
            sage: H2.cusp_normalizing_map(H2.cusps()[1])
            [     2 -1/2*a]
            [     a     -2]
            sage: x=var('x')
            sage: K = NumberField(x^3-36*x-1, names='a')
            sage: H3=HilbertModularGroup(K)
            sage: H3.cusp_normalizing_map(H3.cusps()[1])
            [                     2 1/2*a^2 - 1/2*a - 35/2]
            [                 a + 1                     -8]
            sage: H3.cusp_normalizing_map(H3.cusps()[2])
            [                     3      1/3*a^2 - a - 7/3]
            [1/3*a^2 + 1/3*a - 26/3                      7]
            sage: H3.cusp_normalizing_map(H3.cusps()[3])
            [                     2                      1]
            [1/3*a^2 + 1/3*a - 23/3 1/6*a^2 + 1/6*a - 10/3]

            sage: K4 = NumberField(x**4 - 17*x**2 + 36,names='a'); a=K4.gen()
            sage: H4=HilbertModularGroup(K4)
            sage: H4.cusp_normalizing_map(H4.cusps()[1])
            [                          2 1/4*a^3 - 1/4*a^2 - 4*a + 4]
            [                      a + 1                          -2]

        """
        base_nf = self.base_ring().number_field()
        if not isinstance(cusp, NFCusp) or cusp.number_field() != base_nf:
            raise ValueError(f"Input should be a NF cusp defined over {base_nf}!")
        ca, cb = (cusp.numerator(), cusp.denominator())
        if not (ca, cb) in self._cusp_normalizing_maps:
            # First find the equivalent representative
            # crep, B = self.cusp_representative(cusp,return_map=True)
            # crepa,crepb = crep.numerator(),crep.denominator()
            # crep_normalizing_map = self._cusp_normalizing_maps.get((crepa,crepb))
            # if not crep_normalizing_map:
            # Find a normalizing map of the cusp representative
            a, b, c, d = cusp.ABmatrix()
            det = a * d - b * c
            A = Matrix(self.base_ring().number_field(), 2, 2, [a, b / det, c, d / det])
            # A = B.matrix().inverse()*crep_normalizing_map
            if check:
                infinity = NFCusp(self.base_ring().number_field(), 1, 0)
                if infinity.apply(A.list()) != cusp or A.det() != 1:
                    msg = f"Did not get correct normalizing map A={A} to cusp: {cusp}"
                    raise ArithmeticError(msg)
            logger.debug(f"A={0}".format(A))
            logger.debug("A.det()={0}".format(A.det().complex_embeddings()))
            self._cusp_normalizing_maps_inverse[(ca, cb)] = A.inverse()
            self._cusp_normalizing_maps[(ca, cb)] = A
        if inverse:
            return self._cusp_normalizing_maps_inverse[(ca, cb)]
        else:
            return self._cusp_normalizing_maps[(ca, cb)]

    def apply_cusp_normalizing_map(self, cusp, z, inverse=False):
        """
        Apply the cusp normalising map associated with the cusp to an element z

        INPUT:

        - `cusp` - an instance of NFcusp
        - `z` - an element in
                 - the base number field
                 - the set o cusps
                 -  in ComplexPlaneProductElement__class
        - `inverse` - set to True if applying the inverse map

        EXAMPLES::

            sage: from hilbert_modgroup.all import HilbertModularGroup,ComplexPlaneProductElement
            sage: H2=HilbertModularGroup(10)
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],1.0)
            0.360379610028063
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],1)
            1/6*a - 1/6
            sage: z = NFCusp(H2.base_ring().number_field(),1)
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],z)
            Cusp [a - 1: 6] of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?
            sage: a=H2.base_ring().gens()[1]
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],a)
            3/16*a
            sage: z=ComplexPlaneProductElement([CC(1,0),CC(1,0)]); z
            [1.00000000000000, 1.00000000000000]
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],z)
            [-0.693712943361397, 0.360379610028063]
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],1).complex_embeddings()
            [-0.693712943361397, 0.360379610028063]

            # If we apply a matrix to a an element in K we get back an element in K
            sage: s,t = H2.cusps()[1].numerator(),H2.cusps()[1].denominator()
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],s/t)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: number field element division by zero

            # If we apply the matrix to a cusp we return a cusp.
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],H2.cusps()[1])
            Cusp Infinity of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?

            # Applying the inverse of a cusp normalising map to the same cusp returns infinity.
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],H2.cusps()[1],inverse=True)
            Cusp Infinity of Number Field ... polynomial x^2 - 10 with a = 3.162277660168380?
            sage: c1 = H2.cusps()[1]
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],Infinity)
            1/5*a
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],Infinity) == c1
            False
            sage: q = c1.numerator()/c1.denominator()
            sage: H2.apply_cusp_normalizing_map(H2.cusps()[1],Infinity) == q
            True
            sage: c0, c1 = H2.cusps()
            sage: H2.apply_cusp_normalizing_map(c1,c0) == c1
            True

        """
        a, b, c, d = self.cusp_normalizing_map(cusp, inverse=inverse).list()
        if z == infinity:
            return a / c
        number_field = self.base_ring().number_field()
        if isinstance(z, NFCusp) and z.number_field() == number_field:
            return z.apply([a, b, c, d])
        if z in number_field:
            return (a * z + b) / (c * z + d)
        if isinstance(z, ComplexPlaneProductElement__class) and \
                z.degree() == number_field.absolute_degree():
            return z.apply(matrix(2, 2, [a, b, c, d]))
        raise ValueError("Unsupported type for acting with cusp normalizer! (z={0})".format(z))
