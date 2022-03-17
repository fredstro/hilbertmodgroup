from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.rings.complex_mpc cimport MPComplexField_class, MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber
from sage.categories.map cimport Map

cdef class ComplexPlaneProduct__class(Parent):
    cdef int _degree

cdef class UpperHalfPlaneProduct__class(ComplexPlaneProduct__class):
    pass

cdef class AnytoCPP(Map):
    pass
cdef class ComplexPlaneProductElement__class(Element):
    # cdef double *_x
    # cdef double *_y
    cdef int _degree
    cdef list _x
    cdef list _y
    cdef list _z
    cdef int _prec
    cdef int _verbose
    cdef MPComplexField_class _base_ring

    cdef MPComplexNumber _norm
    cdef RealNumber _imag_norm
    cdef RealNumber _real_norm
    cdef RealNumber _abs_square_norm
    cdef int _norm_set
    cdef int _imag_norm_set
    cdef int _real_norm_set
    cdef int _abs_square_norm_set
    cdef int _is_in_upper_half_plane
    cpdef imag(self)
    cpdef real(self)
    cpdef _is_equal(self, ComplexPlaneProductElement__class other)
    cpdef __copy__(self)
    cpdef z(self)
    cpdef set_prec(self,prec)
    cpdef abs(self)
    cpdef imag_norm(self)
    cpdef real_norm(self)
    cpdef abs_square_norm(self)
    cpdef imag_trace(self)
    cpdef real_trace(self)
    cpdef trace(self)
    cpdef reflect(self)

    cpdef prec(self)
    cpdef norm(self)
    cpdef vector(self)
    cpdef vector_norm(self, p=?)
    cpdef apply(self, m)

cdef class UpperHalfPlaneProductElement__class(ComplexPlaneProductElement__class):
    pass

