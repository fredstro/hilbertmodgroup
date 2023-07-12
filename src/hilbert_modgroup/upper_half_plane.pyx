#cython: language_level=3
#cython: debug=True
r"""

Elements in the upper half plane of degree n

Note: The structure of this class is based on ArithmeticSubgroupElement from sage/modular/arithgroup/arithgroup_element.pyx

"""
import sage
from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
from sage.rings.real_mpfr import RealField
from sage.structure.element cimport Element
from sage.rings.all import Integer, CC
from sage.rings.infinity import Infinity
from sage.structure.parent import Parent
from sage.structure.element cimport parent

from sage.rings.complex_mpfr cimport ComplexNumber
from sage.rings.complex_mpc cimport MPComplexNumber, MPComplexField_class
from sage.rings.complex_mpc import MPComplexField
from cpython.object cimport Py_EQ, Py_NE
from sage.modules.free_module_element import vector


# Constructors for products Complex planes and upper half-planes
def ComplexPlaneProduct(degree, **kwds):
    r"""
    Construct a product of complex planes.

    INPUT:

    - `degree` - integer

    EXAMPLES::

        sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
        sage: ComplexPlaneProduct(2)
        Product of complex planes of degree 2
        sage: ComplexPlaneProduct(3)
        Product of complex planes of degree 3

    """
    if isinstance(degree, sage.rings.number_field.number_field_base.NumberField):
        degree = degree.absolute_degree()
    return ComplexPlaneProduct__class(degree, **kwds)

def UpperHalfPlaneProduct(degree, **kwds):
    r"""
    Construct a product of complex planes.

    INPUT:

    - `degree` - integer

    EXAMPLES::

        sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct
        sage: UpperHalfPlaneProduct(2)
        Product of upper half-planes of degree 2
        sage: UpperHalfPlaneProduct(3)
        Product of upper half-planes of degree 3

    """
    if isinstance(degree, sage.rings.number_field.number_field_base.NumberField):
        degree = degree.absolute_degree()
    return UpperHalfPlaneProduct__class(degree, **kwds)

## Constructors for elements of products of complex planes and upper half-planes
def UpperHalfPlaneProductElement(z, **kwds):
    """
    Construct an element in the product of upper half planes.

    INPUT:

        - ``z`` -- input to construct a tuple of complex number
        - ``kwds`` -- dict.
        - ``degree`` -- positive integer. If a scalar input is given this is the degree of the constructed element.

    OUTPUT:
        - Element of the type UpperHalfPlaneProductElement__class

    EXAMPLES::

        sage: from hilbert_modgroup.all import UpperHalfPlaneProductElement
        sage: UpperHalfPlaneProductElement([1+I,1+I])
        [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
        sage: set(UpperHalfPlaneProductElement([1+I]*10))
        {1.00000000000000 + 1.00000000000000*I}
        sage: len(UpperHalfPlaneProductElement([1+I]*10))
        10
        sage: UpperHalfPlaneProductElement([1,1-I])
        Traceback (most recent call last):
        ...
        ValueError: Point [1.00000000000000, 1.00000000000000 - 1.00000000000000*I] not in upper half-plane!

    """
    if isinstance(z,UpperHalfPlaneProductElement__class):
        parent = kwds.get('parent')
        if parent is z.parent():
            return z
        if parent:
            return UpperHalfPlaneProductElement__class(list(z), parent=parent)
        return z
    prec = kwds.get('prec',getattr(z,'prec',lambda : 53)())
    if hasattr(z,'value'):
        z = z.value()
    if isinstance(z, sage.rings.number_field.number_field_element.NumberFieldElement):
        z = z.complex_embeddings(prec)
    if isinstance(z,list) and not isinstance(z[0], (ComplexNumber, MPComplexNumber)):
       z  = [MPComplexField(prec)(x) for x in z]
    elif not isinstance(z,list) and kwds.get('degree',0)>0:
        z = [MPComplexField(prec)(z)]*kwds.get('degree')
    if 'parent' not in kwds:
        kwds['parent'] = UpperHalfPlaneProduct(degree=len(z))
    return UpperHalfPlaneProductElement__class(z,**kwds)

def ComplexPlaneProductElement(z,**kwds):
    """
    Construct an element in the product of complex planes.

    INPUT:

        - ``z`` -- input to construct a tuple of complex number
        - ``kwds`` -- dict.
            - ``degree`` -- positive integer. If a scalar input is given this is the degree of the constructed element.

    OUTPUT:
        - Element of the type ComplexPlaneProductElement__class

    EXAMPLES::

        sage: from hilbert_modgroup.all import ComplexPlaneProductElement
        sage: z=ComplexPlaneProductElement([CC(1,1),CC(1,1)]); z
        [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
        sage: a=QuadraticField(5).gen()
        sage: ComplexPlaneProductElement(a)
        [-2.23606797749979, 2.23606797749979]
        sage: u0,u1=QuadraticField(5).unit_group().gens()
        sage: ComplexPlaneProductElement(u0)
        [-1.00000000000000, -1.00000000000000]
        sage: ComplexPlaneProductElement(u1)
        [1.61803398874989, -0.618033988749895]

    """
    if isinstance(z,ComplexPlaneProductElement__class):
        parent = kwds.get('parent')
        if parent is z.parent():
            return z
        if parent:
            return ComplexPlaneProductElement__class(list(z), parent=parent)
        return z
    # Get precision in the first hand from kwds, second from z and third set default to 53 bits
    prec = kwds.get('prec',getattr(z,'prec',lambda : 53)())
    if isinstance(z, sage.rings.number_field.number_field_element.NumberFieldElement):
        z = z.complex_embeddings(prec)
    if hasattr(z,'value'):
        z = z.value().complex_embeddings(prec)
    if isinstance(z,list) and not isinstance(z[0],(ComplexNumber,MPComplexNumber)):
        z = [MPComplexField(prec)(x) for x in z]
    elif not isinstance(z,list) and kwds.get('degree',0)>0:
        z = [MPComplexField(prec)(z)]*kwds.get('degree')
    return ComplexPlaneProductElement__class(z,**kwds)


cdef class ComplexPlaneProduct__class(Parent):

    Element = ComplexPlaneProductElement__class

    def __init__(self,degree, **kwds):
        r"""
        Class for a product of complex planes

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct__class
            sage: ComplexPlaneProduct__class(2)
            Product of complex planes of degree 2
            sage: TestSuite(ComplexPlaneProduct__class(2)).run()

        """
        Parent.__init__(self)
        self._degree = degree

    def __hash__(self):
        """
        Return hash of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: hash(ComplexPlaneProduct(2)) == hash('Product of complex planes of degree 2')
            True

        """
        return hash(str(self))

    def construction(self):
        r"""
        No functor exists here but this needs to be defined for coercion to work properly.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: ComplexPlaneProduct(2).construction() is None
            True


        """
        return None

    cpdef _coerce_map_from_(self, S):
        r"""
        Coerce maps from S to self.
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: H=ComplexPlaneProduct(2)
            sage: H._coerce_map_from_(ZZ)
            Generic map:
                From: Integer Ring
                To:   Product of complex planes of degree 2
            sage: H._coerce_map_from_(QuadraticField(5))
            Generic map:
                From: Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?
                To:   Product of complex planes of degree 2
            sage: H._coerce_map_from_(QuadraticField(5))(QuadraticField(5)(1))
             [1.00000000000000, 1.00000000000000]
        """
        if self._coerce_from_hash is None:
            self.init_coerce(False)
        if type(S) == type(self) and S.degree() == self.degree():
            from sage.categories.homset import Hom
            morphism = Hom(self, self).identity()
            morphism._is_coercion = True
            self._coerce_from_hash.set(S, morphism)
            return morphism
        if type(S) == type(self):
            msg = f"Can not coerce UpperHalfPlaneProduct of degree {S.degree()} to degree {self.degree()}"
            raise TypeError(msg)
        try:
            morphism = AnytoCPP(S,self)
            self._coerce_from_hash.set(S, morphism)
            return morphism
        except:
            pass
        return super(ComplexPlaneProduct__class,self)._internal_coerce_map_from(S)

    def _an_element_(self):
        r"""
        Create a typical element of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: ComplexPlaneProduct(2)._an_element_()
            [ - 1.00000000000000*I,  - 1.00000000000000*I]

        """
        return self._element_constructor_([CC(0,-1)]*self.degree())

    def __eq__(self, other):
        r"""
        Check if self is equal to other

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: ComplexPlaneProduct(2) == ComplexPlaneProduct(2)
            True
            sage: ComplexPlaneProduct(2) == ComplexPlaneProduct(3)
            False

        """
        if not isinstance(other,type(self)):
            return False
        return self.degree() == other.degree()

    def __str__(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct__class
            sage: ComplexPlaneProduct__class(2)
            Product of complex planes of degree 2

        """
        return f"Product of complex planes of degree {self._degree}"

    def __repr__(self):
        """
        Representation of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct__class
            sage: ComplexPlaneProduct__class(2)
            Product of complex planes of degree 2

        """
        return str(self)

    def __reduce__(self):
        r"""
        Prepare self for pickling

        TESTS::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: c = ComplexPlaneProduct(2)
            sage: loads(dumps(c)) == c
            True

        """
        return ComplexPlaneProduct, (self.degree(),)

    def degree(self):
        r"""
        Return the degree of this product of complex planes.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProduct
            sage: ComplexPlaneProduct(2).degree()
            2
            sage: ComplexPlaneProduct(3).degree()
            3

        """
        return self._degree

    def _element_constructor_(self,z, **kwds):
        r"""

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProduct
            sage: ComplexPlaneProduct(degree=2)._element_constructor_([1,1])
            [1.00000000000000, 1.00000000000000]
            sage: ComplexPlaneProduct(degree=2)._element_constructor_([1,1+I])
            [1.00000000000000, 1.00000000000000 + 1.00000000000000*I]
        """
        kwds['degree'] = self.degree()
        kwds['parent'] = self
        return ComplexPlaneProductElement(z, **kwds)

    cpdef coerce(self, x):
        r"""
        Coerce x to an element of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProduct
            sage: ComplexPlaneProduct(degree=2).coerce([1,1])
            [1.00000000000000, 1.00000000000000]
            sage: ComplexPlaneProduct(degree=2).coerce([1,1+I])
            [1.00000000000000, 1.00000000000000 + 1.00000000000000*I]

        """
        return self._element_constructor_(x)

from sage.categories.map cimport Map
cdef class AnytoCPP(Map):
    """
    Maps from 'anything' into the class ComplexPlaneProduct

    TODO: implement this as a combination of maps from elements to complex numbers and then to lists.

    """
    cpdef Element _call_(self, x):
        """
        
        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import AnytoCPP, ComplexPlaneProduct
            sage: H = ComplexPlaneProduct(2)
            sage: AnytoCPP(ZZ,H)
            Generic map:
              From: Integer Ring
              To:   Product of complex planes of degree 2
            sage: AnytoCPP(ZZ,H)(1)
            [1.00000000000000, 1.00000000000000]
            sage: AnytoCPP(CC,H)
            Generic map:
              From: Complex Field with 53 bits of precision
              To:   Product of complex planes of degree 2
            sage: AnytoCPP(CC,H)(CC(1,1))
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: AnytoCPP(str,H)("1")
            [1.00000000000000, 1.00000000000000]
            sage: AnytoCPP(str,H)("a")
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'a' to a MPComplexNumber
        """
        cdef ComplexPlaneProduct__class parent = <ComplexPlaneProduct__class>self._codomain
        return parent._element_constructor(x)

    def section(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_mpfr import RRtoRR
            sage: R10 = RealField(10)
            sage: R100 = RealField(100)
            sage: f = RRtoRR(R100, R10)
            sage: f.section()
            Generic map:
              From: Real Field with 10 bits of precision
              To:   Real Field with 100 bits of precision
        """
        return AnytoCPP(self._codomain, self.domain())


cdef class UpperHalfPlaneProduct__class(ComplexPlaneProduct__class):
    r"""
    Class for elements in a product of upper half-planes including the boundary (i.e. imaginary part >=0).

    EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct__class
            sage: UpperHalfPlaneProduct__class(2)
            Product of upper half-planes of degree 2
            sage: TestSuite(UpperHalfPlaneProduct__class(2)).run()

    """
    Element = UpperHalfPlaneProductElement__class

    def _element_constructor_(self,z, **kwds):
        r"""
        Construct an element of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import UpperHalfPlaneProduct
            sage: UpperHalfPlaneProduct(degree=2)._element_constructor_([1,1-I])
            Traceback (most recent call last):
            ...
            ValueError: Point [1.00000000000000, 1.00000000000000 - 1.00000000000000*I] not in upper half-plane!
            sage: UpperHalfPlaneProduct(degree=2)._element_constructor_([1+I,1+2*I])
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 2.00000000000000*I]


        """
        kwds['degree'] = self.degree()
        kwds['parent'] = self
        return UpperHalfPlaneProductElement(z, **kwds)

    def _an_element_(self):
        r"""
        Create a typical element of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct
            sage: UpperHalfPlaneProduct(2)._an_element_()
            [1.00000000000000*I, 1.00000000000000*I]

        """
        return self._element_constructor_([CC(0.0,1.0)]*self.degree())

    def __str__(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct__class
            sage: UpperHalfPlaneProduct__class(2)
            Product of upper half-planes of degree 2

        """
        return f"Product of upper half-planes of degree {self.degree()}"

    def __repr__(self):
        """
        Representation of self.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct__class
            sage: UpperHalfPlaneProduct__class(2)
            Product of upper half-planes of degree 2

        """
        return str(self)

    def __reduce__(self):
        r"""
        Prepare self for pickling

        TESTS::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct
            sage: c = UpperHalfPlaneProduct(2)
            sage: loads(dumps(c)) == c
            True

        """
        return UpperHalfPlaneProduct, (self.degree(),)


cdef class ComplexPlaneProductElement__class(Element):
    r"""
    Class of elements in products of complex planes
    with additional ring structure given by:
    - component-wise multiplication and division
    - multiplication and division by elements of a number field of the same degree as the dimension.

    EXAMPLES::

        sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement__class
        sage: z=ComplexPlaneProductElement__class([CC(1,1),CC(1,1)]); z
        [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
        sage: z.parent()
        Product of complex planes of degree 2
        sage: TestSuite(z).run()
        sage: a=QuadraticField(5).gen()
        sage: ComplexPlaneProductElement__class(a.complex_embeddings())
        [-2.23606797749979, 2.23606797749979]
        sage: u0,u1=QuadraticField(5).unit_group().gens()
        sage: u0 = QuadraticField(5)(u0)
        sage: u1 = QuadraticField(5)(u1)
        sage: ComplexPlaneProductElement__class(u0.complex_embeddings())
        [-1.00000000000000, -1.00000000000000]
        sage: ComplexPlaneProductElement__class(u1.complex_embeddings())
        [1.61803398874989, -0.618033988749895]

    TODO: Inherit from Ring or something? (for speed probably NO!)

    """

    Parent = ComplexPlaneProduct__class

    def __init__(self,zl, verbose=0, *argv, **kwds):
        r"""
        Init self from a list of complex numbers.
        Currently we only work with double (53 bits) precision.

        INPUT:

        - `zl` (list) - list of complex numbers

        EXAMPLES:

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement__class
            sage: from hilbert_modgroup.all import ComplexPlaneProduct,ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement__class([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: H=ComplexPlaneProduct(degree=2)
            sage: ComplexPlaneProductElement([1,2,3],parent=H)
            Traceback (most recent call last):
            ...
            ValueError: Can not construct an element of degree 2 from list of length 3
        """
        self._verbose = verbose
        if verbose>0:
            print("in __init__")
        if not isinstance(zl,list):
            raise ValueError("Need a list to init")
        parent = kwds.get('parent')
        self._degree = len(zl)
        if not parent:
            parent = ComplexPlaneProduct(self._degree)
        if self._degree != parent.degree():
            msg = f"Can not construct an element of degree {parent.degree()} from list of length {len(zl)}"
            raise ValueError(msg)
        if not isinstance(zl[0],(MPComplexNumber,ComplexNumber)):
            raise ValueError("Need a list of MPComplexNumber")
        super().__init__(parent)
        self._prec = zl[0].prec()
        self._base_ring = MPComplexField(self._prec)
        if verbose>0:
            print(zl[0],type(zl[0]))
        if isinstance(zl[0],ComplexNumber):
            self._z = [self._base_ring(z) for z in zl]
        else:
            self._z = zl
        self._x = [z.real() for z in zl]
        self._y = [z.imag() for z in zl]
        if all([x>=0 for x in self._y]):
            self._is_in_upper_half_plane = True
        else:
            self._is_in_upper_half_plane = False
        self._imag_norm = 0.0
        self._real_norm = 0.0
        self._norm = 0.0
        self._abs_square_norm = 0.0
        self._imag_norm_set = 0
        self._real_norm_set = 0
        self._norm_set = 0
        self._abs_square_norm_set = 0


    def _cache_key(self):
        """
        Cache key for self.

        EXAMPLES::

        sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement__class
        sage: z=ComplexPlaneProductElement__class([CC(1,1),CC(2,3)])
        sage: z._cache_key()
        ('ComplexPlaneProductElement__class',
        (1.00000000000000 + 1.00000000000000*I,
        2.00000000000000 + 3.00000000000000*I))

        """
        return (self.__class__.__name__,tuple(self._z))

    def __reduce__(self):
        r"""
        Prepare self for pickling

        TESTS::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement
            sage: c = ComplexPlaneProductElement([1,1])
            sage: loads(dumps(c)) == c
            True

        """
        return ComplexPlaneProductElement, (self.z(),)

    def __hash__(self):
        """
        Hash of self.

        EXAMPLES::

        sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement__class
        sage: z=ComplexPlaneProductElement__class([CC(1,1),CC(2,3)])
        sage: hash(z) # random
        2592654998731023797

        """
        return hash(self._cache_key())

    cpdef z(self):
        r"""
        Return the list of complex numbers in this element.
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement, UpperHalfPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)])
            sage: z.z()
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: z.z()
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z=UpperHalfPlaneProductElement([CC(1,1),CC(2,1)])
            sage: z.z()
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 1.00000000000000*I]

        
        """
        return self._z

    def is_in_upper_half_plane(self):
        """
        Base ring of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: z.is_in_upper_half_plane()
            True
            sage: z=ComplexPlaneProductElement([CC(1,0),CC(2,3)]); z
            [1.00000000000000, 2.00000000000000 + 3.00000000000000*I]
            sage: z.is_in_upper_half_plane()
            True
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.is_in_upper_half_plane()
            False

        """
        return bool(self._is_in_upper_half_plane)

    def as_upper_half_plane_element(self):
        r"""
        Return a copy of self with type UpperHalfPlaneProductElement__class

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: z.as_upper_half_plane_element()
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.as_upper_half_plane_element()
            Traceback (most recent call last):
            ...
            ValueError: Can not convert self to element in product of upper half-planes.
        """
        if not self.is_in_upper_half_plane():
            raise ValueError("Can not convert self to element in product of upper half-planes.")
        return UpperHalfPlaneProductElement(self._z)

    def is_zero(self):
        r"""
        Return true if all components of self is zero

        EXAMPLES:

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(0,0),CC(0,0)])
            sage: z.is_zero()
            True
            sage: z=ComplexPlaneProductElement([CC(0,0),CC(1,0)])
            sage: z.is_zero()
            False
            sage: z=ComplexPlaneProductElement([CC(-1,0),CC(1,0)])
            sage: z.is_zero()
            False
        """

        return all([x==0 for x in self])

    def base_ring(self):
        r"""
        Base ring of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: z.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    cpdef prec(self):
        r"""
        The precision of self.


        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1),CC(4,2)])
            sage: z.prec()
            53
            sage: z=ComplexPlaneProductElement([ComplexField(106)(1,1),ComplexField(106)(2,-1)])
            sage: z.prec()
            106        
        """
        return self._prec

    def degree(self):
        """
        Degree of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: z.degree()
            2

        """
        return self._degree

    def __len__(self):
         """
         Length of self.

         EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: len(z)
            2
         """
         return self.degree()

    cpdef real(self):
        """
        Real parts of self
        
        EXAMPLES::    
    
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z.real()
            [1.00000000000000, 2.00000000000000]
            
        """
        return self._x

    cpdef imag(self):
        """
        Imaginary parts of self


        EXAMPLES::    
    
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z.imag()
            [1.00000000000000, 3.00000000000000]
        """
        return self._y

    cpdef __copy__(self):
        """
        Copy self.
        
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: w=copy(z); w
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: w==z
            True
        """
        return self.__class__(self._z, verbose=self._verbose)

    def __repr__(self):
        """
        String representation of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]

        """
        return str(list(self))

    def __getitem__(self,i):
        """
        Get the items of self.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z[0]
            1.00000000000000 + 1.00000000000000*I
            sage: z[1]
            2.00000000000000 - 1.00000000000000*I
        """
        if isinstance(i, (int, Integer)) and 0 <= i < self._degree:
            return self._z[i]
        else:
            raise IndexError

    cpdef _is_equal(self, ComplexPlaneProductElement__class other):
        """
        Return 1 if ``self`` is equal to ``other``
        
        INPUT:
        - ``other`` - Element of the type ``ComplexPlaneProductElement__class``
        
        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: w1=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: w2=ComplexPlaneProductElement([CC(1,1),CC(2,1)])
            sage: z._is_equal(z)
            1
            sage: z._is_equal(w1)
            1
            sage: z._is_equal(w2)
            0

        """
        cdef int i
        for i in range(self.degree()):
            if self._x[i] != other._x[i] or self._y[i] != other._y[i]:
                return 0
        return 1

    def __richcmp__(self, right, int op):
        """
        Compare self with other

        INPUT:
        - `right`
        - `op`

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: w=copy(z); w
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: w==z
            True
            sage: w!=z
            False
            sage: w>z
            Traceback (most recent call last):
            ...
            NotImplementedError: Ordering of points in H^n is not implemented!
        """
        res=1
        if op != Py_EQ and op != Py_NE:
            raise NotImplementedError("Ordering of points in H^n is not implemented!")
        if type(self) != type(right) or right.degree() != self.degree():
            res=0
        else:
            res = self._is_equal(right)
        if op == Py_NE:
            res = 1 - res
        return bool(res)

    cpdef imag_norm(self):
        """
        Return the product of all imaginary parts of self.
        

        EXAMPLES::
            
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.imag_norm()
            -1.00000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(0,0)]); z
            [1.00000000000000 + 1.00000000000000*I, 0]
            sage: z.imag_norm()
            0.000000000000000
        """
        cdef int i
        if self._imag_norm_set==0:
            self._imag_norm=self.base_ring().base_ring()(1)
            for i in range(self._degree):
                self._imag_norm = self._imag_norm*self._y[i]
            self._imag_norm_set=1
        return self._imag_norm

    cpdef abs(self):
        """
        Return the element consisting of the absolute value of all elements.


        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.abs() # abs tol 1e-10
            [1.41421356237310, 2.23606797749979]
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(0,0)]); z
            [1.00000000000000 + 1.00000000000000*I, 0]
            sage: z.imag_norm()
            0.000000000000000
        """
        return ComplexPlaneProductElement([abs(self._z[i]) for i in range(self._degree)])

    cpdef abs_square_norm(self):
        r"""
        Return the norm of |z|^2 
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: z.abs_square_norm() # abs tol 1e-10
            10.0000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(0,0)])
            sage: z.abs_square_norm() # abs tol 1e-10
            0.000000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(1,0)])
            sage: z.abs_square_norm() # abs tol 1e-10
            2.000000000000000
            
        """
        cdef int i
        if not self._abs_square_norm_set:
            self._abs_square_norm = self.base_ring().base_ring()(1)
            for i in range(self._degree):
                self._abs_square_norm = self._abs_square_norm * (self._x[i]**2 + self._y[i]**2)
            self._abs_square_norm_set = 1
        return self._abs_square_norm

    cpdef real_norm(self):
        """
        Return the product of all real parts of self.
        
        EXAMPLES::
            
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.real_norm()
            2.00000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(0,0)]); z
            [1.00000000000000 + 1.00000000000000*I, 0]
            sage: z.real_norm()
            0.000000000000000


        """
        cdef int i
        if not self._real_norm_set:
            self._real_norm=self.base_ring().base_ring()(1)
            for i in range(self._degree):
                self._real_norm = self._real_norm*self._x[i]
            self._real_norm_set=1
        return self._real_norm

    cpdef norm(self):
        """
        Return the product of all components of self.
        
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.norm()
            3.00000000000000 + 1.00000000000000*I
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(0,0)]); z
            [1.00000000000000 + 1.00000000000000*I, 0]
            sage: z.norm()
            0
            sage: z=ComplexPlaneProductElement([ComplexField(106)(1,1),ComplexField(106)(2,-1)]); z
            [1.000000000000000000000000000000 + 1.000000000000000000000000000000*I, 2.000000000000000000000000000000 - 1.000000000000000000000000000000*I]
            sage: z.norm()
            3.000000000000000000000000000000 + 1.000000000000000000000000000000*I
        """
        cdef int i
        if self._norm_set==0:
            self._norm = self.base_ring()(1)
            for i in range(self._degree):
                self._norm = self._norm*self._z[i]
            self._norm_set=1
        return self._norm

    cpdef vector(self):
        r"""
        Return self as a vector
        
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.vector()
            (1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I)
            sage: z=ComplexPlaneProductElement([ComplexField(106)(1,1),ComplexField(106)(2,-1)]); z
            [1.000000000000000000000000000000 + 1.000000000000000000000000000000*I, 2.000000000000000000000000000000 - 1.000000000000000000000000000000*I]
            sage: z.vector()
            (1.000000000000000000000000000000 + 1.000000000000000000000000000000*I, 2.000000000000000000000000000000 - 1.000000000000000000000000000000*I)
            sage: type(z.vector())
            <class 'sage.modules.free_module_element.FreeModuleElement_generic_dense'>
            sage: z.vector().base_ring()
            Complex Field with 106 bits of precision

        """
        return vector(self)

    cpdef vector_norm(self,p=2):
        r"""
        Return the Euclidean norm of self as a vector in C^n
        
        INPUT:
        - `p` (integer) default = 2 (L2-norm). Other options include =1 (L1-norm) or =0 (Infinity-norm)
        
        Note: This is about twice as fast as doing z.vector().norm()

        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.vector_norm()
            2.64575131106459
            sage: z.vector_norm(Infinity)
            2.23606797749979
            sage: z.vector_norm(1)
            3.65028153987288
            sage: z.vector_norm() == z.vector().norm()
            True
            sage: z.vector_norm(5)
            2.27959492535969
        """
        from sage.rings.infinity import infinity
        cdef int i
        if p == infinity:
            return max([abs(z) for z in self._z])
        p = self.base_ring().base_ring()(p)
        res = sum([abs(z)**p for z in self._z])
        if p != 1:
            res = res**(p**-1)
        return res

    def change_ring(self,R):
        """
        Change the base ring of self.

        INPUT:
        - `R` -- MOComplexField

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.change_ring(MPComplexField(106));z
            [1.000000000000000000000000000000 + 1.000000000000000000000000000000*I, 2.000000000000000000000000000000 - 1.000000000000000000000000000000*I]
            sage: z.change_ring(MPComplexField(22));z
            [1.00000 + 1.00000*I, 2.00000 - 1.00000*I]

        """
        if not isinstance(R,MPComplexField_class):
            raise ValueError(f"Can not coerce self into {R}")
        self.set_prec(R.prec())

    cpdef set_prec(self,prec):
        """
        Change the precision of self.
        
        INPUT:
        - `prec`
        
        EXAMPLES::
        
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.set_prec(106);z
            [1.000000000000000000000000000000 + 1.000000000000000000000000000000*I, 2.000000000000000000000000000000 - 1.000000000000000000000000000000*I]
            sage: z.set_prec(22);z
            [1.00000 + 1.00000*I, 2.00000 - 1.00000*I]

        """
        self._base_ring = MPComplexField(prec)
        self._z = [self._base_ring(z) for z in self]
        self._x = [self._base_ring.base_ring()(x) for x in self._x]
        self._y = [self._base_ring.base_ring()(y) for y in self._y]

    @classmethod
    def _extract_left_right_parent(cls,left,right):
        """
        Convert the argument left and right to elements of the type ComplexPlaneProduct_class

        INPUT:

        - ``left`` -- object
        - ``right`` -- object

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement,ComplexPlaneProductElement__class
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: w=ComplexPlaneProductElement([CC(1,1),CC(2,1)])
            sage: ComplexPlaneProductElement__class._extract_left_right_parent(z,w)
            ([1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I],
             [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 1.00000000000000*I],
             Product of complex planes of degree 2)

        """
        if isinstance(left, ComplexPlaneProductElement__class):
            right = ComplexPlaneProductElement(right, prec=left.prec(), degree=left.degree())
            parent = left.parent()
        elif isinstance(right, ComplexPlaneProductElement__class):
            left = ComplexPlaneProductElement(left, prec=right.prec(), degree=right.degree())
            parent = right.parent()
        else:
            raise ValueError("One of left or right must be of the type ComplexPlaneProductElement__class!")
        return left,right,parent

    cdef _add_(self, other):
        """
        Add ``other`` to ``self`` and convert to ``parent``. Used by the ``__add__`` method.
        
        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement,ComplexPlaneProductElement__class     
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement__class       
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: w=ComplexPlaneProductElement([CC(2,2),CC(3,3)]); w
            [2.00000000000000 + 2.00000000000000*I, 3.00000000000000 + 3.00000000000000*I]
            sage: z+w
            [3.00000000000000 + 3.00000000000000*I, 5.00000000000000 + 6.00000000000000*I]
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z+1
            [2.00000000000000 + 1.00000000000000*I, 3.00000000000000 - 1.00000000000000*I]
            sage: K=QuadraticField(5)
            sage: a=K.gen()
            sage: z+a
            [-1.23606797749979 + 1.00000000000000*I, 4.23606797749979 - 1.00000000000000*I]
            sage: a+z
            [-1.23606797749979 + 1.00000000000000*I, 4.23606797749979 - 1.00000000000000*I]
            sage: w=ComplexPlaneProductElement([CC(2,2),CC(3,3),CC(4,4)]); w
            [2.00000000000000 + 2.00000000000000*I, 3.00000000000000 + 3.00000000000000*I, 4.00000000000000 + 4.00000000000000*I]
            sage: z+w
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +: 'Product of complex planes of degree 2' and 'Product of complex planes of degree 3'
            sage: from hilbert_modgroup.all import UpperHalfPlaneProduct
            sage: z=UpperHalfPlaneProduct(degree=2)([1+I,1+2*I]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 2.00000000000000*I]
            sage: z+2
            [3.00000000000000 + 1.00000000000000*I, 3.00000000000000 + 2.00000000000000*I]
            sage: 2+z
            [3.00000000000000 + 1.00000000000000*I, 3.00000000000000 + 2.00000000000000*I]
            sage: type(z+2) == type(z)
            True
            sage: type(2+z) == type(z)
            True

        """
        if self._degree != other.degree() or self._prec != other.prec():
            raise TypeError
        return self._parent([self._z[i] + other[i] for i in range(self.degree())])

    cdef _neg_(self):
        """
        Negative of self.
        
        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)])
            sage: -z
            [-1.00000000000000 - 1.00000000000000*I, -2.00000000000000 - 3.00000000000000*I]
       
        """
        return self._parent([-self._z[i] for i in range(self.degree())])

    cdef _sub_(self,other):
        """
        Subtract ``other`` from ``self``

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,3)]); z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 3.00000000000000*I]
            sage: w=ComplexPlaneProductElement([CC(2,2),CC(3,3)]); w
            [2.00000000000000 + 2.00000000000000*I, 3.00000000000000 + 3.00000000000000*I]
            sage: z+w
            [3.00000000000000 + 3.00000000000000*I, 5.00000000000000 + 6.00000000000000*I]
            sage: z-w
            [-1.00000000000000 - 1.00000000000000*I, -1.00000000000000]
            sage: z-1
            [1.00000000000000*I, 1.00000000000000 + 3.00000000000000*I]
            sage: K=QuadraticField(5)
            sage: a=K.gen()
            sage: z-a
            [3.23606797749979 + 1.00000000000000*I, -0.236067977499790 + 3.00000000000000*I]
            sage: a-z
            [-3.23606797749979 - 1.00000000000000*I, 0.236067977499790 - 3.00000000000000*I]
            sage: w=ComplexPlaneProductElement([CC(2,2),CC(3,3),CC(4,4)]); w
            [2.00000000000000 + 2.00000000000000*I, 3.00000000000000 + 3.00000000000000*I, 4.00000000000000 + 4.00000000000000*I]
            sage: z-w
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for -: 'Product of complex planes of degree 2' and 'Product of complex planes of degree 3'


        """
        if self._degree != other.degree() or self._prec != other.prec():
            raise TypeError
        # Try to make an element of the same class as self and if it doesn't work, coerce to complex plane product element
        try:
            return self._parent([self._z[i] - other[i] for i in range(self.degree())])
        except ValueError:
            return ComplexPlaneProductElement__class([self._z[i] - other[i] for i in range(self.degree())])

    cdef _mul_(self, other):
        r"""
        Multiply ``self`` by ``other`` and convert to ``parent``.
        
        INPUT:
        - `other` - element of product of complex planes
        - `parent` - a parent class 
        
        EXAMPLES::
            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement,ComplexPlaneProductElement__class
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: w=ComplexPlaneProductElement([CC(2,1),CC(1,1)])
            sage: z*w
            [1.00000000000000 + 3.00000000000000*I, 3.00000000000000 + 1.00000000000000*I]
            sage: w*z
            [1.00000000000000 + 3.00000000000000*I, 3.00000000000000 + 1.00000000000000*I]
            sage: z*z
            [2.00000000000000*I, 3.00000000000000 - 4.00000000000000*I]
            sage: w*w
            [3.00000000000000 + 4.00000000000000*I, 2.00000000000000*I]
            sage: K=QuadraticField(5)
            sage: a=K.gen()
            sage: z*a
            [-2.23606797749979 - 2.23606797749979*I, 4.47213595499958 - 2.23606797749979*I]
            sage: a*z
            [-2.23606797749979 - 2.23606797749979*I, 4.47213595499958 - 2.23606797749979*I]
            sage: z*CC(1,3)
            [-2.00000000000000 + 4.00000000000000*I, 5.00000000000000 + 5.00000000000000*I]
            sage: z*5
            [5.00000000000000 + 5.00000000000000*I, 10.0000000000000 - 5.00000000000000*I]
            sage: u0,u1=K.unit_group().gens()
            sage: z*u1
            [1.61803398874989 + 1.61803398874989*I, -1.23606797749979 + 0.618033988749895*I]

            # check Upper half plane elements
            sage: from hilbert_modgroup.all import UpperHalfPlaneProduct
            sage: z=UpperHalfPlaneProduct(degree=2)([1+I,1+2*I]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 2.00000000000000*I]
            sage: z*2
            [2.00000000000000 + 2.00000000000000*I, 2.00000000000000 + 4.00000000000000*I]
            sage: 2*z
            [2.00000000000000 + 2.00000000000000*I, 2.00000000000000 + 4.00000000000000*I]
            sage: type(z*2) == type(z)
            True
            sage: type(2*z) == type(z)
            True
            sage: -1*z
            [-1.00000000000000 - 1.00000000000000*I, -1.00000000000000 - 2.00000000000000*I]
            sage: w=UpperHalfPlaneProduct(degree=2)([-1+I,2+I]); w
            [-1.00000000000000 + 1.00000000000000*I, 2.00000000000000 + 1.00000000000000*I]
            sage: w*z
            [-2.00000000000000, 5.00000000000000*I]

        """
        if self._degree != other.degree() or self._prec != other.prec():
            raise TypeError
        try:
            new_element = [self._z[i]*other[i] for i in range(self.degree())]
            return self._parent(new_element)
        except ValueError:
            return ComplexPlaneProductElement__class(new_element)

    cdef _div_(self, other):
        r"""
        Divide self by other.
        
        INPUT:
        - `other` - element of product of complex planes
        - `parent` - parent class
        
        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import ComplexPlaneProductElement,ComplexPlaneProductElement__class
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: w=ComplexPlaneProductElement([CC(2,1),CC(1,1)])
            sage: z / w
            [0.600000000000000 + 0.200000000000000*I, 0.500000000000000 - 1.50000000000000*I]
            sage: w / z
            [1.50000000000000 - 0.500000000000000*I, 0.200000000000000 + 0.600000000000000*I]
            sage: z / z
            [1.00000000000000, 1.00000000000000]
            sage: w / (z / z)
            [2.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: w=ComplexPlaneProductElement([CC(2,1),CC(1,1)])
            sage: K=QuadraticField(5)
            sage: a=K.gen()
            sage: z/a
            [-0.447213595499958 - 0.447213595499958*I, 0.894427190999916 - 0.447213595499958*I]
            sage: a/z
            [-1.11803398874989 + 1.11803398874989*I, 0.894427190999916 + 0.447213595499958*I]
            sage: z/2
            [0.500000000000000 + 0.500000000000000*I, 1.00000000000000 - 0.500000000000000*I]
            sage: z/[0,1]
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Can not divide by zero!

            # check Upper half plane elements
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProduct
            sage: z=UpperHalfPlaneProduct(degree=2)([1+I,1+2*I]); z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 2.00000000000000*I]
            sage: z/2
            [0.500000000000000 + 0.500000000000000*I, 0.500000000000000 + 1.00000000000000*I]
            sage: type(z/2) == type(z)
            True
            sage: -1/z
            [-0.500000000000000 + 0.500000000000000*I, -0.200000000000000 + 0.400000000000000*I]
            sage: type(-1/z)==type(z)
            True
            sage: 1/z
            [0.500000000000000 - 0.500000000000000*I, 0.200000000000000 - 0.400000000000000*I]

        """
        if self._degree != other.degree() or self._prec != other.prec():
            raise TypeError
        if any([z == 0 for z in other]):
            raise ZeroDivisionError("Can not divide by zero!")
        new_element = [self._z[i]/other[i] for i in range(self.degree())]
        try:
            return self.parent()(new_element)
        except ValueError:
            return ComplexPlaneProductElement__class(new_element)


    def __pow__(self, power, modulo):
        """
        Self to the power 'power' defined component-wise:
        If `power` is a scalar:
                z^power = (z_1^power,...,z_n^power)
        If `power` is an element of this class:
                z^power = (z_1^power_1,...,z_n^power_n)

        INPUT:

        - `power` -- complex number (will be coerced to the base_ring of self)
        - `modulo` -- dummy argument (ignored)

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)])
            sage: z**2
            [2.00000000000000*I, 3.00000000000000 - 4.00000000000000*I]
            sage: z**[2,2]
            [2.00000000000000*I, 3.00000000000000 - 4.00000000000000*I]
            sage: z**[1,2]
            [1.00000000000000 + 1.00000000000000*I, 3.00000000000000 - 4.00000000000000*I]
            sage: w=ComplexPlaneProductElement([CC(2,1),CC(1,1)])
            sage: z**w
            [-0.309743504928494 + 0.857658012588736*I, 3.35025931507288 + 1.18915022150040*I]
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(0,0)]);
            sage: z**2
            [2.00000000000000*I, 0]
            sage: z**-2
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Can not divide component by 0!
            sage: z**[-2,1]
            [ - 0.500000000000000*I, 0]
            sage: z**[-2,0]
            [ - 0.500000000000000*I, 1.00000000000000]
            """

        if modulo:
            raise RuntimeError("__pow__ dummy argument ignored")
        try:
            if not isinstance(power,(list,ComplexPlaneProductElement__class)):
                power = [power]*self.degree()
            power = self.parent()(power)
            if any(power[i].real() < 0 and self[i] == 0 for i in range(self.degree())):
                raise ZeroDivisionError("Can not divide component by 0!")
            new_element = [z**power[i] for i,z in enumerate(self)]
            try:
                return self.parent()(new_element)
            except ValueError:
                return ComplexPlaneProductElement__class(new_element)
        except TypeError as e:
            raise TypeError(f"Can not coerce {power} to self.base_ring(). {e}")


    cpdef trace(self):
        """
        Trace of self, i.e. sum over all components.


        EXAMPLES::
                
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]);z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.trace()
            3.00000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1),CC(4,2)])
            sage: z.trace()
            7.00000000000000 + 2.00000000000000*I
        
        """
        return sum(self)

    cpdef real_trace(self):
        """
        Trace of real part self, i.e. sum over all components.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]);z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.real_trace()
            3.00000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1),CC(4,2)])
            sage: z.real_trace()
            7.00000000000000
        
        
        """
        return sum(self.real())

    cpdef imag_trace(self):
        """
        Trace of imaginary part of self, i.e. sum over all components.

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1)]);z
            [1.00000000000000 + 1.00000000000000*I, 2.00000000000000 - 1.00000000000000*I]
            sage: z.imag_trace()
            0.000000000000000
            sage: z=ComplexPlaneProductElement([CC(1,1),CC(2,-1),CC(4,2)])
            sage: z.imag_trace()
            2.00000000000000

        """
        return sum(self.imag())

    cpdef apply(self, m):
        r"""
        Apply the matrix m to self. 
    
        INPUT:
        
        - ``m`` -- matrix
        
        EXAMPLES::
            
            sage: from hilbert_modgroup.all import ComplexPlaneProductElement, UpperHalfPlaneProductElement, HilbertModularGroup
            sage: z=ComplexPlaneProductElement([CC(3,1),CC(1,-1)]);z
            [3.00000000000000 + 1.00000000000000*I, 1.00000000000000 - 1.00000000000000*I]
            sage: z=ComplexPlaneProductElement([CC(3,1),CC(-1,1)]);z
            [3.00000000000000 + 1.00000000000000*I, -1.00000000000000 + 1.00000000000000*I]
            sage: N = matrix(ZZ,[[1,0],[0,1]])
            sage: z.apply(N)
            [3.00000000000000 + 1.00000000000000*I, -1.00000000000000 + 1.00000000000000*I]        
            sage: A=matrix(ZZ,[[0,-1],[1,0]])
            sage: z.apply(A)
            [-0.300000000000000 + 0.100000000000000*I, 0.500000000000000 + 0.500000000000000*I]
            sage: H5 = HilbertModularGroup(5)
            sage: A = H5(A)
            sage: z=UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement__class
            sage: z.apply(A)
            [-0.300000000000000 + 0.100000000000000*I, 0.500000000000000 + 0.500000000000000*I]            
            sage: isinstance(_, UpperHalfPlaneProductElement__class)
            True
            sage: a=H5.base_ring().number_field().gen()
            sage: A = H5.cusp_normalizing_map(NFCusp(H5.base_ring().number_field(),a,1+a))
            sage: z.apply(A)
            [1.67855780465319 + 0.0271280367431345*I, 0.717919314224174 + 0.0871677256896697*I]
    
        """
        try:
            aa, bb, cc, dd = m.list()
        except (AttributeError,ValueError):
            raise ValueError("Need a 2 x 2 matrix or object that contains a list of 4 elements to act on self.")

        try:
            aa, bb, cc, dd = m.list()
        except (AttributeError, ValueError):
            raise ValueError(
                "Need a 2 x 2 matrix or object that contains a list of 4 elements to act on self.")
        if hasattr(aa,'complex_embeddings'):
            a = aa.complex_embeddings()
            b = bb.complex_embeddings()
            c = cc.complex_embeddings()
            d = dd.complex_embeddings()
        else:
            a = [RealField()(aa)]*self._degree
            b = [RealField()(bb)]*self._degree
            c = [RealField()(cc)]*self._degree
            d = [RealField()(dd)]*self._degree
        ## Component - wise application of map
        denominators = [(c[i]*z+d[i]) for i, z in enumerate(self._z)]
        if 0 in denominators:
            return Infinity
        zlist = [ (a[i]*z+b[i])/denominators[i] for i, z in enumerate(self._z)]
        return self.parent()(zlist)

    def as_ComplexPlaneProductElement(self):
        """
        Convert self to an element in the product of complex planes.

        EXAMPLES::

            sage: from hilbert_modgroup.all import UpperHalfPlaneProductElement
            sage: z = UpperHalfPlaneProductElement([1+I,1+I])
            sage: z.as_ComplexPlaneProductElement()
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]

            sage: z = UpperHalfPlaneProductElement([1+I,1+I,1+I,1+I])
            sage: z.as_ComplexPlaneProductElement()
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]

        """
        if self.__class__ == ComplexPlaneProductElement__class:
            return self
        return ComplexPlaneProductElement(list(self))

    def permute(self,s):
        r"""
        The permutation group acts by permuting the components.

        EXAMPLES::

            sage: from hilbert_modgroup.all import UpperHalfPlaneProductElement
            sage: from sage.groups.perm_gps.constructor import PermutationGroupElement
            sage: g=PermutationGroupElement([2,1])
            sage: z = UpperHalfPlaneProductElement([1+I,2+2*I])
            sage: z.permute(g)
            [2.00000000000000 + 2.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: z = UpperHalfPlaneProductElement([1+I,2+I,3+I,4+I])
            sage: g=PermutationGroupElement([2,1])
            sage: z.permute(g)
            [2.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I, 3.00000000000000 + 1.00000000000000*I, 4.00000000000000 + 1.00000000000000*I]
            sage: g=PermutationGroupElement([2,1,4,3])
            sage: z.permute(g)
            [2.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I, 4.00000000000000 + 1.00000000000000*I, 3.00000000000000 + 1.00000000000000*I]


        """
        if not is_PermutationGroupElement(s):
            raise ValueError("Input must be a permutation group element")
        znew = [0 for i in range(self._degree)]
        for i in range(self._degree):
            si = s(i+1)-1
            znew[si]=self.z()[i]
        return self.parent()(znew)

    cpdef reflect(self):
        r"""
        Reflection of self in each coordinate in the imaginary axis.

        EXAMPLES::

            sage: from hilbert_modgroup.all import UpperHalfPlaneProductElement
            sage: z = UpperHalfPlaneProductElement([1+I,2+2*I])
            sage: z.reflect()
            [-1.00000000000000 + 1.00000000000000*I, -2.00000000000000 + 2.00000000000000*I]

            sage: z = UpperHalfPlaneProductElement([1+I,1+I,1+I,1+I])
            sage: z.reflect()
            [-1.00000000000000 + 1.00000000000000*I, -1.00000000000000 + 1.00000000000000*I, -1.00000000000000 + 1.00000000000000*I, -1.00000000000000 + 1.00000000000000*I]

        """
        znew = [0 for i in range(self._degree)]
        for i in range(self.degree()):
            znew[i] = MPComplexField(self._prec)(-self._x[i],self._y[i])
        return self.parent()(znew)



cdef class UpperHalfPlaneProductElement__class(ComplexPlaneProductElement__class):
    r"""
    Class for elements of products of the upper half plane.

    """

    def __init__(self, zl, verbose=0, *argv, **kwds):
        r"""
        Init self from list of complex numbers.

        EXAMPLES::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement__class
            sage: from sage.rings.complex_mpc import MPComplexField
            sage: MPC = MPComplexField(53)
            sage: z = UpperHalfPlaneProductElement__class([MPC(1,1),MPC(1,1)])
            sage: TestSuite(z).run()
            sage: z
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]
            sage: UpperHalfPlaneProductElement__class([MPC(1,1),MPC(1,1),MPC(1,1)])
            [1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I, 1.00000000000000 + 1.00000000000000*I]

        TESTS:
            sage: UpperHalfPlaneProductElement__class([1+I,1+I])
            Traceback (most recent call last):
            ...
            ValueError: Need a list of MPComplexNumber

        """
        super().__init__(zl,verbose=verbose,*argv,**kwds)
        if not self.is_in_upper_half_plane():
            raise ValueError("Point {0} not in upper half-plane!".format(zl))

    def __reduce__(self):
        r"""
        Prepare self for pickling

        TESTS::

            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
            sage: c = UpperHalfPlaneProductElement([1,1])
            sage: loads(dumps(c)) == c
            True

        """
        return UpperHalfPlaneProductElement, (self.z(),)

    def imag_log(self):
        r"""
        The logarithmic embedding of the products of upper half-planes into R^k given by log.

        EXAMPLES:

            sage: from hilbert_modgroup.all import UpperHalfPlaneProductElement
            sage: z=UpperHalfPlaneProductElement([CC(1,3),CC(1,2)])
            sage: z.imag_log()
            (1.09861228866811, 0.693147180559945)

        """
        return vector([x.imag().log() for x in self])

    def hyp_dist(self, w, dtype=1):
        r"""
        Return "a" hyperbolic distance between ``self`` and ``w``

        INPUT:
        - `w` -- element of the upper half plane
        - `dtype` -- integer (default=1).
                     If dtype == 0 the combined distance is the sum of all components
                     If dtype == 1 the combined distance is the maximum of all components.

        EXAMPLES::

            sage: from hilbert_modgroup.all import UpperHalfPlaneProductElement
            sage: z = UpperHalfPlaneProductElement([1+I,2+2*I])
            sage: w = UpperHalfPlaneProductElement([2+I,2+2*I])
            sage: z.hyp_dist(z)
            0
            sage: z.hyp_dist(w)
            0.962423650119207

        """
        if dtype == 1:
            maxd = 0
        if not isinstance(w,UpperHalfPlaneProductElement__class) or w.degree() != self.degree():
            raise ValueError("w must be an element of the same degree")

        if w.z() == self.z():
            return 0
        distances = []
        for i in range(self.degree()):
            ab1 = abs(self.z()[i] - w.z()[i])
            ab2 = abs(self.z()[i] - MPComplexField(self._prec)(w.real()[i], -w.imag()[i]))
            l = ((ab1 + ab2) / (ab2 - ab1)).log()
            distances.append(l)

        if dtype == 0:
            return sum(distances)
        else:
            return max(distances)

    cpdef apply(self, m):
        r"""
        Apply the matrix m to self. 

        INPUT:

        - ``m`` -- matrix

        EXAMPLES::

            sage: from hilbert_modgroup.all import ComplexPlaneProductElement, UpperHalfPlaneProductElement, HilbertModularGroup
            sage: H5 = HilbertModularGroup(5)
            sage: A=matrix(ZZ,[[0,-1],[1,0]])
            sage: A = H5(A)
            sage: z=UpperHalfPlaneProductElement([CC(3,1),CC(-1,1)])
            sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement__class
            sage: z.apply(A)
            [-0.300000000000000 + 0.100000000000000*I, 0.500000000000000 + 0.500000000000000*I]            
            sage: isinstance(_, UpperHalfPlaneProductElement__class)
            True
            sage: a=H5.base_ring().number_field().gen()
            sage: A = H5.cusp_normalizing_map(NFCusp(H5.base_ring().number_field(),a,1+a))
            sage: z.apply(A)
            [1.67855780465319 + 0.0271280367431345*I, 0.717919314224174 + 0.0871677256896697*I]

        """
        new_point = super(UpperHalfPlaneProductElement__class,self).apply(m)
        return new_point.as_upper_half_plane_element()