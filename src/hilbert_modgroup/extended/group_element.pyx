from sage.structure.element cimport MultiplicativeGroupElement
from sage.structure.richcmp cimport richcmp
from sage.rings.all import ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.misc.cachefunc import cached_method
from hilbert_modgroup.upper_half_plane cimport ComplexPlaneProductElement__class

cdef class ExtendedHilbertModularGroupElement(MultiplicativeGroupElement):
    cdef Matrix_generic_dense __x

    def __init__(self, parent, x, check=True):
        """
        Create an element of a extended hilbert Modular Group.

        INPUT:

        - ``parent`` -- an arithmetic subgroup

        - `x` -- data defining a 2x2 matrix
                 which lives in parent

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: TestSuite(H).run()
            sage: x, y = H.OK().gens()
            sage: ExtendedHilbertModularGroupElement(H, [1, x, 0, 1])
            [          1 1/2*a + 1/2]
            [          0           1]
            sage: ExtendedHilbertModularGroupElement(H,[1,x,0,1]) in H
            True
            sage: K3.<a> = QuadraticField(3)
            sage: lattice_ideal = K3.different()
            sage: f = K3.unit_group().fundamental_units()[0]
            sage: H = ExtendedHilbertModularGroup(K3, tp_units = False)
            sage: A = ExtendedHilbertModularGroupElement(H, [1, 0, 0, f]); A
            Traceback (most recent call last):
            ...
            TypeError: matrix testing must have determinant 1
            sage: H = ExtendedHilbertModularGroup(K3, tp_units = True)
            sage: A = ExtendedHilbertModularGroupElement(H, [1, 0, 0, -f]); A
            Traceback (most recent call last):
            ...
            TypeError: matrix must have determinant equal to totally positive unit
        """
        from .group_class import ExtendedHilbertModularGroup_class
        if not isinstance(parent, ExtendedHilbertModularGroup_class):
            raise ValueError("parent (= {0}) must be a Extended Hilbert Modular group".format(parent))
        x = MatrixSpace(parent.base_ring(), 2, 2)(x, copy=True, coerce=True)
        if parent.tp_units():
            if not (x.determinant().is_unit() and x.determinant().is_totally_positive()):
                raise TypeError("matrix must have determinant equal to totally positive unit")
        else:
            if not (x.determinant() == 1):
                raise TypeError("matrix testing must have determinant 1")
        if not x[0, 1] in parent.lattice_ideal().inverse():
            raise TypeError("matrix must have M[0, 1] entry in lattice_ideal.inverse()")
        if not x[1, 0] in (parent.lattice_ideal() * parent.level_ideal()):
            raise TypeError("matrix must have M[1, 0] element in lattice_ideal*level_ideal")
        if not (x[0, 0] in parent.OK() and x[1, 1] in parent.OK()):
            raise TypeError("matrix must have elements M[0, 0], M[1, 1] in OK")
        x.set_immutable()
        MultiplicativeGroupElement.__init__(self, parent)
        self.__x = x

    def __iter__(self):
        r"""
        Iterate over the entries of ``self``.

        The entries are yielded in row-major order: upper left, upper right,
        lower left, lower right.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: x, y = H.OK().gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, x, 0, 1])
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
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: level_ideal = K2.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.list()
            [1, 1/10*a + 1/2, 0, 1]
        """
        return list(self)

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::
         sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
         sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
         sage: K2.<a> = QuadraticField(5)
         sage: lattice_ideal = K2.different()
         sage: level_ideal = K2.fractional_ideal(3)
         sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal); H
         Hilbert modular group PGL_2^+(Fractional ideal (1) + Fractional ideal (-a), Fractional ideal (3)) of
         order 2 over Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
         consisting of matrices of determinant in U_K^+.

        """
        return repr(self.__x)

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: level_ideal = K2.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: latex(A)
            \left(\begin{array}{rr}
             1 & \frac{1}{10} \sqrt{5} + \frac{1}{2} \\
             0 & 1
            \end{array}\right)
        """
        return self.__x._latex_()

    cpdef _richcmp_(self, right_r, int op):
        r"""
        Compare ``self`` to ``right_r``, where ``right_r`` is guaranteed to
        have the same parent as ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K3.<a> = QuadraticField(3)
            sage: H = ExtendedHilbertModularGroup(K3)
            sage: x, y = H.OK().gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, x, 0, 1])
            sage: B = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A == B
            False
            sage: A != B
            True
            sage: A == A
            True

        """
        cdef ExtendedHilbertModularGroupElement right = <ExtendedHilbertModularGroupElement> right_r
        return richcmp(self.__x, right.__x, op)

    def __nonzero__(self):
        r"""
        Return ``True``, since ``self`` is an element of a group of
        invertible matrices and is never the zero matrix.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: x, y = H.OK().gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, x, 0, 1])
            sage: bool(A)
            True

        """
        return True

    def __invert__(self):
        """
        Return the inverse of matrix associate to self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: level_ideal = K2.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.__invert__()
            [            1 -1/10*a - 1/2]
            [            0             1]
        """
        det = self.__x.determinant()
        return self.parent([self.__x[1, 1] / det, -self.__x[0, 1] / det,
                            -self.__x[1, 0] / det, self.__x[0, 0] / det])

    def matrix(self):
        """
        Return the matrix of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: level_ideal = K2.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.matrix()
            [           1 1/10*a + 1/2]
            [           0            1]
        """
        return self.__x

    def number_field(self):
        """
        Return the number field of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: level_ideal = K2.fractional_ideal(3)
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal, level_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.number_field()
            Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
        """
        return self.parent().number_field()

    def determinant(self):
        """
        Return the determinant of ``self``, which is always a positive unit.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.determinant()
            1
        """
        return self.matrix().determinant()

    def a(self):
        """
        Return the upper left entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.a()
            1

        """
        return self.__x[0, 0]

    def b(self):
        """
        Return the upper right entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.b()
            1/10*a + 1/2

        """
        return self.__x[0, 1]

    def c(self):
        """
        Return the lower left entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.c()
            0

        """
        return self.__x[1, 0]

    def d(self):
        """
        Return the lower right entry of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: A.d()
            1

        """
        return self.__x[1, 1]

    cpdef _mul_(self, right):
        """
        Return the self*right

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K3.<a> = QuadraticField(3)
            sage: H = ExtendedHilbertModularGroup(K3)
            sage: x, y = H.OK().gens()
            sage: A = ExtendedHilbertModularGroupElement(H, [1, x, 0, 1])
            sage: B = ExtendedHilbertModularGroupElement(H, [1, y, 0, 1])
            sage: C = A*B
            sage: C
            [    1 a + 1]
            [    0     1]
        """
        return self.__class__(self.parent(), self.__x * (<ExtendedHilbertModularGroupElement> right).__x, check=False)

    def acton(self, z):
        """
        Return the result of the action of ``self`` on z as a fractional linear
        transformation.

        EXAMPLES::
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2, lattice_ideal)
            sage: x, y = (lattice_ideal**-1).gens()
            sage: f = K2.unit_group().fundamental_units()[0]
            sage: A = ExtendedHilbertModularGroupElement(H, [1, 0, y**-1, f**2 ]); A
            [           1            0]
            [-1/2*a + 5/2 -1/2*a + 3/2]


        An example of A acting on a symbolic variable::

            sage: z = var('z')
            sage: A.acton(z)
            -2*z/(z*(sqrt(5) - 5) + sqrt(5) - 3)

        An example of A acting on an element of the base field:

            sage: a = K2.gen()
            sage: z = a/7 +1/12; z
            1/7*a + 1/12
            sage: z in A.number_field()
            True
            sage: A.acton(z)
            2941/23362*a + 3449/23362

        An example with complex numbers::

            sage: C.<i> = ComplexField()
            sage: A.acton(i)
            0.672251363384123 + 0.185805707042680*I

        An example with the cusp infinity::

            sage: A.acton(infinity)
            1/10*a + 1/2

        An example acting on the NFCusp elements
            sage: from hilbert_modgroup.extended.all import NFCusp_wrt_lattice_ideal
            sage: c = NFCusp_wrt_lattice_ideal(lattice_ideal, -a, 1); c
            Cusp [-a: 1] of Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790? with
            respect to  lattice_ideal
            sage: A.acton(c)
            4/29*a + 15/29

        """
        from sage.rings.infinity import InfinityElement, infinity
        if isinstance(z, InfinityElement):
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
        except Exception:
            raise ValueError(f"Can not apply self to z of type: {type(z)}")

    cpdef _acton_complex_plane_element(self, ComplexPlaneProductElement__class z):
        r"""
        Act on an element of type :class:`ComplexPlaneProductElement__class`
        via the Moebius transformation defined by the complex embeddings
        of ``self``.

        INPUT:

        - ``z`` -- an instance of :class:`ComplexPlaneProductElement__class`

        OUTPUT:

        An element of the same type as ``z``, obtained by applying
        the Moebius transformation componentwise using the complex
        embeddings of ``self``.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: S = ExtendedHilbertModularGroupElement(H, [0, -1, 1, 0])
            sage: z = UpperHalfPlaneProductElement([CC(0,1), CC(0,1)]); z
            [1.00000000000000*I, 1.00000000000000*I]
            sage: S._acton_complex_plane_element(z)
            [1.00000000000000*I, 1.00000000000000*I]

        """
        result = []
        if len(z) != len(self.complex_embeddings()):
            raise ValueError("Need element of the same degree!")
        CF = z.base_ring()
        for i, Aemb in enumerate(self.complex_embeddings()):
            a, b, c, d = [CF(x) for x in Aemb.list()]
            result.append((a * z[i] + b) / (c * z[i] + d))
        return z.parent()(result)

    def __getitem__(self, q):
        r"""
        Fetch entries by direct indexing.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: A = ExtendedHilbertModularGroupElement(H, [3, 2, 1, 1])
            sage: A[0, 0]
            3
        """
        return self.__x[q]

    def __hash__(self):
        r"""
        return the hash value associated to self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: A = ExtendedHilbertModularGroupElement(H, [3, 2, 1, 1]); A
            [3 2]
            [1 1]
            sage: hash(A)
            4303834509739506107
        """
        return hash(self.__x)

    def __reduce__(self):
        """
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: f = K2.unit_group().fundamental_units()[0]
            sage: A = ExtendedHilbertModularGroupElement(H, [1, 0, 0, f**2]); A
            [           1            0]
            [           0 -1/2*a + 3/2]
            sage: A.__reduce__()
            (Hilbert modular group PGL_2^+(Fractional ideal (1) + Fractional ideal (1), Fractional ideal (1))
            of order 2 over Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
            consisting of matrices of determinant in U_K^+.,
            (
            [           1            0]
            [           0 -1/2*a + 3/2]
            ))
        """
        return self.parent(), (self.__x,)

    def trace(self):
        r"""
        Return the trace of the trace of the matrix of self.

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: K2.<a> = QuadraticField(5)
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: f = K2.unit_group().fundamental_units()[0]
            sage: A = ExtendedHilbertModularGroupElement(H, [1, 0, 0, f**2]); A
            [           1            0]
            [           0 -1/2*a + 3/2]
            sage: A.matrix().trace()
            -1/2*a + 5/2
            sage: A.trace()
            5

        """
        return self.matrix().trace().trace()

    @cached_method
    def complex_embeddings(self, prec=53):
        """
        Return a list of matrices which are entry-wise complex embeddings of self

        INPUT:
        - ``prec`` integer (default=53) number of bits precision

        EXAMPLES::

            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroupElement
            sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup
            sage: K2.<a> = QuadraticField(5)
            sage: lattice_ideal = K2.different()
            sage: H = ExtendedHilbertModularGroup(K2)
            sage: f = K2.unit_group().fundamental_units()[0]
            sage: A = ExtendedHilbertModularGroupElement(H, [1, 0, 0, f**2]); A
            [           1            0]
            [           0 -1/2*a + 3/2]
            sage: A.complex_embeddings()
            [
            [ 1.00000000000000 0.000000000000000]
            [0.000000000000000  2.61803398874989],
            <BLANKLINE>
            [ 1.00000000000000 0.000000000000000]
            [0.000000000000000 0.381966011250105]
            ]
            sage: A.complex_embeddings(103)
            [
            [ 1.00000000000000000000000000000 0.000000000000000000000000000000]
            [0.000000000000000000000000000000  2.61803398874989484820458683437],
            <BLANKLINE>
            [ 1.00000000000000000000000000000 0.000000000000000000000000000000]
            [0.000000000000000000000000000000 0.381966011250105151795413165634]
            ]
        """
        emb_a = self.a().complex_embeddings(prec)
        emb_b = self.b().complex_embeddings(prec)
        emb_c = self.c().complex_embeddings(prec)
        emb_d = self.d().complex_embeddings(prec)
        M = MatrixSpace(emb_a[0].parent(), 2, 2)
        return [M([emb_a[i], emb_b[i], emb_c[i], emb_d[i]]) for i in range(len(emb_a))]
