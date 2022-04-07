"""
Element of Hilbert modular groups


AUTHORS:

- Fredrik Stromberg (2021)

NOTE: The structure of this class is based on ArithmeticSubgroupElement from sage/modular/arithgroup/arithgroup_element.pyx
"""
from sage.structure.element cimport MultiplicativeGroupElement
from sage.structure.richcmp cimport richcmp
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.misc.cachefunc import cached_method

from hilbert_modgroup.upper_half_plane cimport ComplexPlaneProductElement__class, UpperHalfPlaneProductElement__class

cdef class HilbertModularGroupElement(MultiplicativeGroupElement):


    cdef Matrix_generic_dense __x

    def __init__(self, parent, x, check=True):
        """
        Create an element of a hilbert Modular Group.

        INPUT:

        - ``parent`` -- an arithmetic subgroup

        - `x` -- data defining a 2x2 matrix over ZZ
                 which lives in parent


        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: from hilbert_modgroup.hilbert_modular_group_element import HilbertModularGroupElement
            sage: H=HilbertModularGroup(5)
            sage: TestSuite(H).run()
            sage: x,y=H.base_ring().gens()
            sage: HilbertModularGroupElement(H,[1,x,0,1])
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: HilbertModularGroupElement(H,[1,y,0,1])
            [1 a]
            [0 1]
            sage: HilbertModularGroupElement(H,[1,y,0,1]) in H
            True

        """
        if not 'HilbertModularGroup_class' in parent.__class__.__name__:
            raise TypeError("parent (= {0}) must be a Hilbert Modular group".format(parent))
        x = MatrixSpace(parent.base_ring(),2,2)(x, copy=True, coerce=True)
        if x.determinant() != 1:
            raise TypeError("matrix must have determinant 1")
        x.set_immutable()
        MultiplicativeGroupElement.__init__(self, parent)
        self.__x = x

    def __iter__(self):
        """
        Iterate over self.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: list(A)
            [1, 1/2*a + 1/2, 0, 1]

        """
        yield self.__x[0, 0]
        yield self.__x[0, 1]
        yield self.__x[1, 0]
        yield self.__x[1, 1]

    def list(self):
        r"""
        List of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: A.list()
            [1, 1/2*a + 1/2, 0, 1]
        """
        return list(self)

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

        sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
        sage: H=HilbertModularGroup(5); H
        Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?

        """
        return repr(self.__x)

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: latex(A)
            \left(\begin{array}{rr}
            1 & \frac{1}{2} \sqrt{5} + \frac{1}{2} \\
            0 & 1
            \end{array}\right)

        """
        return self.__x._latex_()

    cpdef _richcmp_(self, right_r, int op):
        r"""
        Compare self to right, where right is guaranteed to have the same
        parent as self.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: B = H([1,y,0,1])
            sage: A == B
            False
            sage: A != B
            True
            sage: A < B
            True
            sage: A > B
            False

        """
        cdef HilbertModularGroupElement right = <HilbertModularGroupElement>right_r
        return richcmp(self.__x, right.__x, op)

    def __nonzero__(self):
        """
        Return ``True``, since the ``self`` lives in SL(2,\Z), which does not
        contain the zero matrix.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: bool(A)
            True
        """
        return True

    cpdef _mul_(self, right):
        """
        Return self * right.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: B = H([1,y,0,1])
            sage: C = A*B
            sage: C
            [          1 3/2*a + 1/2]
            [          0           1]
            sage: C.parent()
            Hilbert Modular Group ... x^2 - 5 with a = 2.236067977499790?
          
        """
        return self.__class__(self.parent(), self.__x * (<HilbertModularGroupElement> right).__x, check=False)

    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: H([1,x,0,1]).__invert__()
            [           1 -1/2*a - 1/2]
            [           0            1]
        """
        return self._parent(
                [self.__x.get_unsafe(1, 1), -self.__x.get_unsafe(0, 1),
                 -self.__x.get_unsafe(1, 0), self.__x.get_unsafe(0, 0)]
                )

    def matrix(self):
        """
        Return the matrix corresponding to ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1]); A
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: A.matrix()
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: type(A.matrix())
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        """
        return self.__x

    def determinant(self):
        """
        Return the determinant of ``self``, which is always 1.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1])
            sage: A.determinant()
            1
        """
        return self.base_ring().one()

    def det(self):
        """
        Return the determinant of ``self``, which is always 1.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A =H([1,x,0,1])
            sage: A.determinant()
            1
        """
        return self.determinant()

    def a(self):
        """
        Return the upper left entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1]); A
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: A.a()
            1
        """
        return self.__x.get_unsafe(0, 0)

    def b(self):
        """
        Return the upper right entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1]); A
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: A.b()
            1/2*a + 1/2
        """
        return self.__x.get_unsafe(0, 1)

    def c(self):
        """
        Return the lower left entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1]); A
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: A.c()
            0

        """
        return self.__x.get_unsafe(1, 0)

    def d(self):
        """
        Return the lower right entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,0,1]); A
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: A.d()
            1

        """
        return self.__x.get_unsafe(1, 1)

    def acton(self, z):
        """
        Return the result of the action of ``self`` on z as a fractional linear
        transformation.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,x,x+2]); A
            [          1 1/2*a + 1/2]
            [1/2*a + 1/2 1/2*a + 5/2]


        An example of A acting on a symbolic variable::

            sage: z = var('z')
            sage: A.acton(z)
            (2*z + sqrt(5) + 1)/(z*(sqrt(5) + 1) + sqrt(5) + 5)

        An example of A acting on an element of the base field:

            sage: K = A.base_ring().number_field(); K
            Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
            sage: a = K.gen()
            sage: z = a/7 +1/12; z
            1/7*a + 1/12
            sage: z in A.base_ring().number_field()
            True
            sage: A.acton(z)
            21251/92062*a - 3947/92062

        An example with complex numbers::

            sage: C.<i> = ComplexField()
            sage: A.acton(i)
            0.475683661041614 + 0.0636610018750175*I

        An example with the cusp infinity::

            sage: A.acton(infinity)
            1/2*a - 1/2

        An example which maps a finite cusp to infinity::

            sage: A.acton(-a)
            +Infinity

        An example acting on the NFCusp elements

            sage: c=NFCusp(K,-a,1)
            sage: c
            Cusp [-a: 1] of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
            sage: A.acton(c)
            +Infinity

        another example of acting on a NFCusp element

            sage: c=NFCusp(K,a/7 +1/12); c
            Cusp [12*a + 7: 84] of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
            sage: A.acton(c)
            21251/92062*a - 3947/92062

        Example acting on points in the upper half-plane

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)]); z
            [1.00000000000000*I, 1.00000000000000*I]
            sage: A.acton(z)
            [-0.642350327708281 + 0.436338998124982*I, 0.475683661041614 + 0.0636610018750175*I]
            sage: S=H([0,-1,1,0])
            sage: S.acton(z)
            [1.00000000000000*I, 1.00000000000000*I]

        NOTE: when acting on instances of cusps the return value
        is still an element of the underlying number field or infinity (Note the presence of
        '+', which does not show up for cusp instances)::


        TESTS:

        We cover the remaining case, i.e., infinity mapped to infinity::

            sage: H([1, K.gen(), 0, 1]).acton(infinity)
            +Infinity
        """
        from sage.rings.infinity import is_Infinite, infinity
        if is_Infinite(z):
            if self.c() != 0:
                return self.a() / self.c()
            else:
                return infinity
        if hasattr(z, 'denominator') and hasattr(z, 'numerator'):
            p = z.numerator()
            q = z.denominator()
            P = self.a() * p + self.b() * q
            Q = self.c() * p + self.d() * q
            if not Q and P:
                return infinity
            else:
                return P / Q
        if isinstance(z, ComplexPlaneProductElement__class):
            return self._acton_complex_plane_element(z)
        try:
            return (self.a() * z + self.b()) / (self.c() * z + self.d())
        except:
            raise ValueError(f"Can not apply self to z of type: {type(z)}")

    cpdef _acton_complex_plane_element(self, ComplexPlaneProductElement__class z):
        """
        Act on an element of the type ComplexPlaneProductElement__class
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: from hilbert_modgroup.all import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,x,x+2]); A
            [          1 1/2*a + 1/2]
            [1/2*a + 1/2 1/2*a + 5/2]
            sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)]); z
            [1.00000000000000*I, 1.00000000000000*I]
            sage: A._acton_complex_plane_element(z)
            [-0.642350327708281 + 0.436338998124982*I, 0.475683661041614 + 0.0636610018750175*I]
            sage: S=H([0,-1,1,0])
            sage: S._acton_complex_plane_element(z)
            [1.00000000000000*I, 1.00000000000000*I]

        
        """
        result = []
        if len(z) != len(self.complex_embeddings()):
            raise ValueError("Need element of the same degree!")
        for i, Aemb in enumerate(self.complex_embeddings()):
            a, b, c, d = Aemb.list()
            result.append((a*z[i] + b)/(c*z[i]+d))
        return z.parent()(result)

    def __getitem__(self, q):
        r"""
        Fetch entries by direct indexing.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: H([3,2,1,1])[0,0]
            3
        """
        return self.__x[q]

    def __hash__(self):
        r"""
        Return a hash value.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,x,x+2]); hash(A)
            -3779633680309204004
        """
        return hash(self.__x)

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLES::

            sage: (SL2Z.1).__reduce__()
            (Modular Group SL(2,Z), (
            [1 1]
            [0 1]
            ))
        """
        return self.parent(), (self.__x,)

    def trace(self):
        r"""
        Return the trace of the trace of the matrix of self.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([1,x,x,x+2]); A
            [          1 1/2*a + 1/2]
            [1/2*a + 1/2 1/2*a + 5/2]
            sage: A.matrix().trace()
            1/2*a + 7/2
            sage: A.trace()
            7

        """
        return self.matrix().trace().trace()

    def multiplicative_order(self, check=True):
        r"""
        Return the multiplicative order of this element.

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: H.one().multiplicative_order()
            1
            sage: H([-1,0,0,-1]).multiplicative_order()
            2
            sage: H([0,-1,1,0]).multiplicative_order()
            4
            sage: H([0,-1,1,1]).multiplicative_order()
            6
            sage: H([1,1,0,1]).multiplicative_order()
            +Infinity
            sage: x,y=H.base_ring().gens()
            sage: A = H([0,x**-1,-x,0]); A
            [           0  1/2*a - 1/2]
            [-1/2*a - 1/2            0]
            sage: A.multiplicative_order()
            4
            sage: A=H([0,-1,1,x]); A
            [          0          -1]
            [          1 1/2*a + 1/2]
            sage: A.multiplicative_order()
            10
            sage: A = H([1,x,x,x+2]); A
            [          1 1/2*a + 1/2]
            [1/2*a + 1/2 1/2*a + 5/2]
            sage: A.multiplicative_order()
            +Infinity

        """
        m = self.matrix()
        if m.is_one():
            return ZZ(1)
        elif (-m).is_one():
            return ZZ(2)
        t = m.trace()
        # We check if t is in the list of possible orders for this number field
        order = self.parent().order_of_elliptic_element_of_trace(t)
        if check:
            max_order = order or max(self.parent().orders_of_elliptic_elements())
            for i in range(2,max_order+1):
                m = m*self.matrix()
                if m.is_one() and (not order or (i < order or not m.is_one() and i == order)):
                    raise ArithmeticError("Could not find order of element!")
        if not order:
            return infinity
        return order

    @cached_method
    def complex_embeddings(self, prec=53):
        """
        Return a list of matrices which are entry-wise complex embeddings of self

        INPUT:
        - ``prec`` integer (default=53) number of bits precision

        EXAMPLES::

            sage: from hilbert_modgroup.hilbert_modular_group_class import HilbertModularGroup
            sage: H=HilbertModularGroup(5)
            sage: x,y=H.base_ring().gens()
            sage: A = H([0,x**-1,-x,0]); A
            [           0  1/2*a - 1/2]
            [-1/2*a - 1/2            0]
            sage: A.complex_embeddings()
            [
            [0.000000000000000 -1.61803398874989]
            [0.618033988749895 0.000000000000000],
            [0.000000000000000 0.618033988749895]
            [-1.61803398874989 0.000000000000000]
            ]
            sage: A.complex_embeddings(103)
            [
            [0.000000000000000000000000000000 -1.61803398874989484820458683437]
            [0.618033988749894848204586834366 0.000000000000000000000000000000],
            [0.000000000000000000000000000000 0.618033988749894848204586834366]
            [-1.61803398874989484820458683437 0.000000000000000000000000000000]
            ]
        """
        emb_a = self.a().complex_embeddings(prec)
        emb_b = self.b().complex_embeddings(prec)
        emb_c = self.c().complex_embeddings(prec)
        emb_d = self.d().complex_embeddings(prec)
        M = MatrixSpace(emb_a[0].parent(),2,2)
        return [ M([emb_a[i], emb_b[i], emb_c[i], emb_d[i]]) for i in range(len(emb_a))]

