#cython: profile=True
"""
Cython versions of Pullback algorithms.

Note: These are the algorithms that needs optimizing to make it all faster.
"""
from sage.all import ceil,floor,cartesian_product
from sage.modular.cusps_nf import NFCusp
from sage.structure.sage_object cimport SageObject
from hilbert_modgroup.upper_half_plane cimport UpperHalfPlaneProductElement__class
# from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from sage.rings.number_field.number_field_element cimport NumberFieldElement
from sage.rings.number_field.number_field_element import is_NumberFieldElement


cpdef integral_coordinates_in_box(bounds):
    """
    Find all integral points inside a box in R^n.

    INPUT:
    - `bounds` -- a list of upper and lower bounds defining a parallelepiped of positive volume.

    EXAMPLES::
    
        sage: from hilbert_modgroup.pullback_cython import integral_coordinates_in_box
        sage: integral_coordinates_in_box([(-1,1),(-1,1)])
        [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]        
        sage: integral_coordinates_in_box([(-1,1),(-1.5,1.5)])
        [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]        
        sage: integral_coordinates_in_box([(-1,1),(-1,1),(-0.5,0.5)])
        [(-1, -1, 0),
         (-1, 0, 0),
         (-1, 1, 0),
         (0, -1, 0),
         (0, 0, 0),
         (0, 1, 0),
         (1, -1, 0),
         (1, 0, 0),
         (1, 1, 0)]
        sage: integral_coordinates_in_box([(-1,1),(-1,1),(-0.5,0.5)])
        [(-1, -1, 0),
         (-1, 0, 0),
         (-1, 1, 0),
         (0, -1, 0),
         (0, 0, 0),
         (0, 1, 0),
         (1, -1, 0),
         (1, 0, 0),
         (1, 1, 0)]
    """
    if not isinstance(bounds, list):
        raise ValueError
    n = len(bounds)
    coordinates = []
    for i in range(n):
        lower = ceil(bounds[i][0])
        upper = floor(bounds[i][1])
        if lower > upper:
            raise ValueError("Bounds must give interval of positive length.")
        coordinates.append(range(lower, upper + 1))
    # result = [tuple(x) for x in cartesian_product(coordinates)]
    # return result
    return list(cartesian_product(coordinates))

cpdef lattice_elements_in_box(lattice_basis, lattice_bounds, coordinate_bounds):
    """
    Return all coordinates of elements in a lattice in R^n restricted to a specific box and with
    coordinates in another box.

    EXAMPLES::
    
        sage: from hilbert_modgroup.pullback_cython import lattice_elements_in_box
        sage: lattice_elements_in_box([[0,1],[1,0]],[(-1,1),(-1,1)],[(-1,1),(-1,1)])
        [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
        sage: lattice_elements_in_box([[2,0],[0,2]],[(-2,2),(-1,1)],[(-1,1),(-1,1)])
        [(-1, 0), (0, 0), (1, 0)]
        sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
        sage: P1 = HilbertPullback(HilbertModularGroup(5))
        sage: lb, ib = P1._get_lattice_and_ideal_basis()
        sage: lattice_elements_in_box(lb,[(-2,2),(-2,2)],[(-1,1),(-1,1)])
        [(-1, -1), (-1, 0), (0, -1), (0, 0), (0, 1), (1, 0), (1, 1)]
        sage: lattice_elements_in_box(lb,[(-2,2),(-2,2)],[(-1,1),(0,1)])
        [(-1, 0), (0, 0), (0, 1), (1, 0), (1, 1)]
        
    """
    coordinates = integral_coordinates_in_box(coordinate_bounds)
    result = []
    n = len(lattice_basis[0])
    for coordinate_vector in coordinates:
        is_within_bounds = True
        for i in range(n):
            alpha_i = 0.0
            for j in range(n):
                alpha_i = alpha_i + lattice_basis[i][j] * coordinate_vector[j]
            if alpha_i < lattice_bounds[i][0] or alpha_i > lattice_bounds[i][1]:
                # We need to discard this
                is_within_bounds = False
                break
        # If we are within the bounds we add the number field element.
        if is_within_bounds:
            result.append(coordinate_vector)
    return result

cpdef coordinates_to_ideal_elements(coordinates,ideal_basis):
    r"""
    Return elements of an ideal given by coordinates and a basis.
    
    INPUT:
    - `coordinates` -- a list of coordinates (as tuples or lists) of the same length as the basis.
    - `ideal_basis` -- the basis of an ideal.
    
    EXAMPLES::
    
        sage: from hilbert_modgroup.pullback_cython import coordinates_to_ideal_elements
        sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback
        sage: P1 = HilbertPullback(HilbertModularGroup(5))
        sage: lb, ib = P1._get_lattice_and_ideal_basis()
        sage: coordinates_to_ideal_elements([[1,1],(2,3)],ib)
        [1/2*a + 1/2, 3/2*a + 1/2]
        
    TESTS::
        
        sage: coordinates_to_ideal_elements([[1,1],(2,2,1)],ib)
        Traceback (most recent call last):
        ...
        ValueError: Coordinate need to have same length as basis!
                        
    """
    result = []
    n = len(ideal_basis)
    for coordinate_vector in coordinates:
        if len(coordinate_vector) != n:
            raise ValueError("Coordinate need to have same length as basis!")
        element = 0
        for i,b in enumerate(ideal_basis):
            element += b * coordinate_vector[i]
        result.append(element)
    return result


cpdef find_candidate_cusps(p, z):
    r"""
    Return candidates for closest cusp to the point ``z`` in the upper half-plane. 
    
    INPUT:
    - `p` -- object of type HilbertPullback
    - `z` -- point in the upper half-plane
    
    OUTPUT:
    - list of tuples (sigma,rho) representing cusps 
    
    EXAMPLES::
    
        sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
        sage: from hilbert_modgroup.pullback_cython import find_candidate_cusps
        sage: H1 = HilbertModularGroup(5)
        sage: P1 = HilbertPullback(H1)
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: find_candidate_cusps(P1,z)
        [(1, 0), (-1, -1), (1, -1), (0, 1)]
        sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
        sage: find_candidate_cusps(P1,z)
        [(1, 0), (-1, -1), (1, -1), (-1, -1/2*a + 1/2), (1, -1/2*a + 1/2), (0, 1)]
        
        sage: H2 = HilbertModularGroup(10)
        sage: P2 = HilbertPullback(H2)
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: find_candidate_cusps(P2,z)
         [(1, 0),
         (-4, -4),
         (-3, -4),
         ...
         (3, -a),
         (4, -a),
         (0, 1)]
         sage: len(_)
         100
        sage: z=UpperHalfPlaneProductElement([CC(2.58,0.5),CC(0.5,0.5)]) # long time (1 second)
        sage: len(find_candidate_cusps(P2,z)) # long time (1 second)
        271
        
    """
    ideal_basis = p._construct_ideal(1).integral_basis()
    lattice_basis = p.basis_matrix_ideal()
    n = len(lattice_basis[0])
    # Make lattice basis to a nested list to avoid creation of FreeModule elements
    lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
    coordinate_bounds = p._bound_for_sigma_coordinates(z)
    embedding_bounds = p._bound_for_sigma_embeddings(z)
    coordinate_bounds = [(-b,b) for b in coordinate_bounds]
    embedding_bounds = [(-b,b) for b in embedding_bounds]
    sigma_candidates_coordinates = lattice_elements_in_box(lattice_basis,
                                               embedding_bounds,
                                               coordinate_bounds)
    sigma_candidates = coordinates_to_ideal_elements(sigma_candidates_coordinates,
                                                             ideal_basis)
    result = [(p.number_field()(1),p.number_field()(0))]
    for s in sigma_candidates:
        if s == 0:
            continue
        rho_coordinate_bounds = p._bound_for_rho_coordinates(z, s)
        rho_coordinate_bounds = [(-b,b) for b in rho_coordinate_bounds]
        rho_embedding_bounds = p._bound_for_rho_embeddings(z, s)

        rho_candidates_coordinates = lattice_elements_in_box(lattice_basis,
                                                 rho_embedding_bounds,
                                                 rho_coordinate_bounds)
        rho_candidates = coordinates_to_ideal_elements(rho_candidates_coordinates,
                                                         ideal_basis)

        for r in rho_candidates:
            if (r,s) in result:
                continue
            if (-r,-s) in result:
                continue
            if r == 0 and ((0,1) in result or s != 1):
               continue
            result.append((r,s))
    return result

cpdef find_closest_cusp(p, z, return_multiple=False):
    r"""
    Return the closest cusp(s) to z
    
    INPUT:
    - ``p`` -- A HilbertPullback object
    - ``z`` -- point in the upper half-plane
    - ``return_multiple`` -- boolean (default False), set to True to return multiple closest cusp if it is no unique. 
    
    OUTPUT:
    - tuple representing a cusp (if return_multiple == False)
    - list of tuples representing cusps (if return_multiple == True) 
    
    EXAMPLES::

        sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
        sage: from hilbert_modgroup.pullback_cython import find_closest_cusp
        sage: H1 = HilbertModularGroup(5)
        sage: P1 = HilbertPullback(H1)
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: find_closest_cusp(P1,z)
        (1, 0) 
        sage: find_closest_cusp(P1,z,return_multiple=True)
        [(1, 0), (0, 1)]
        sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
        sage: find_closest_cusp(P1,z)
        (0, 1)
        sage: find_closest_cusp(P1,z,return_multiple=True)
        [(0, 1)]
        
        sage: H2 = HilbertModularGroup(10)
        sage: P2 = HilbertPullback(H2)
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: find_closest_cusp(P2,z)
        (1, 0)
        sage: find_closest_cusp(P2,z,return_multiple=True)
        [(1, 0), (0, 1)]
        sage: z=UpperHalfPlaneProductElement([CC(-2.58,0.5),CC(0.5,0.5)])
        sage: find_closest_cusp(P2,z,return_multiple=True)
        [(a - 5, -a)]
        
    """
    cusp_candidates = find_candidate_cusps(p, z)
    min_cusp = cusp_candidates[0]
    min_d = distance_to_cusp(p,min_cusp[0], min_cusp[1], z)
    if return_multiple:
        min_cusp = [min_cusp]
    for c in cusp_candidates[1:]:
        d = distance_to_cusp(p, c[0],c[1], z)
        if d < min_d:
            if return_multiple:
                min_cusp = [c]
            else:
                min_cusp = c
        if d == min_d and return_multiple:
            min_cusp.append(c)
    return min_cusp

cpdef distance_to_cusp(SageObject p, NumberFieldElement r, NumberFieldElement s,
                        UpperHalfPlaneProductElement__class z):
    r"""
    Return the distance from an element in the upper half-plane to a cusp. 
    
    INPUT:
    - ``p`` - Should be an object of type HilbertPullback (but this is not recognized by Cython so I use SageObject instead)
    - ``r`` - Number field element, numerator of cusp
    - ``s`` - Number field element, denomnator of cusp
    - ``z`` - Point in the upper half-plane
    
    EXAMPLES::
    
        sage: from hilbert_modgroup.all import HilbertModularGroup, HilbertPullback, UpperHalfPlaneProductElement
        sage: from hilbert_modgroup.pullback_cython import distance_to_cusp
        sage: H1 = HilbertModularGroup(5)
        sage: P1 = HilbertPullback(H1)
        sage: K = P1.number_field()
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: distance_to_cusp(P1,K(0),K(1),z) # abs tol 1e-10
        1.0
        sage: distance_to_cusp(P1,K(1),K(1),z) # abs tol 1e-10
        2.0
        sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
        sage: distance_to_cusp(P1,K(0),K(1),z) # abs tol 1e-10
        0.707106781186548        
        sage: H2 = HilbertModularGroup(10)
        sage: P2 = HilbertPullback(H2)
        sage: K = P2.number_field()
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: distance_to_cusp(P2,K(0),K(1),z) # abs tol 1e-10
        1.0
        sage: z=UpperHalfPlaneProductElement([CC(-2.58,0.5),CC(0.5,0.5)])
        sage: a = K.gen()
        sage: distance_to_cusp(P2,a-5,-a,z) # abs tol 1e-10
        1.01308408502928
 
    """
    cdef list rlist,slist
    ideal_rs = r.parent().ideal(r,s)
    rlist = r.complex_embeddings()
    slist = s.complex_embeddings()
    n = len(slist)
    d = 1
    for i in range(n):
        d = d* ((-z._x[i] * slist[i] + rlist[i]) ** 2 * z._y[i] ** -1 + slist[i] ** 2 * z._y[i])
    d = ideal_rs.norm()**-1*d.sqrt()
    return d