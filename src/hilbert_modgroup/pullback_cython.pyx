#cython: profile=True
"""
Cython versions of Pullback algorithms.

Note: These are the algorithms that needs optimizing to make it all faster.
"""
from sage.all import ceil,floor,cartesian_product
from sage.arith.misc import gcd
from sage.modular.cusps_nf import NFCusp
from sage.rings.infinity import Infinity
from sage.structure.sage_object cimport SageObject
from hilbert_modgroup.upper_half_plane cimport UpperHalfPlaneProductElement__class
# from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from sage.rings.number_field.number_field_element cimport NumberFieldElement
from sage.rings.number_field.number_field_element import is_NumberFieldElement
from sage.matrix.constructor import identity_matrix
from sage.matrix.all import Matrix

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

cpdef lattice_elements_in_box(lattice_basis, lattice_bounds, coordinate_bounds, norm_bound=None):
    """
    Return all coordinates of elements in a lattice in R^n restricted to a specific box and with
    coordinates in another box.

    INPUT:
    - ``lattice_basis`` -- a list of lists of a basis of the embedded lattice
    - ``lattice_bounds`` -- the bounds for the embeddings
    - ``coordinate_bounds`` -- the bounds for the coordinates
    - ``norm_bound`` -- bounds for the norm (default None - meaning that no bounds are applied)
    
    EXAMPLES::
    
        sage: from hilbert_modgroup.pullback_cython import lattice_elements_in_box
        sage: lattice_elements_in_box([[0,1],[1,0]],[(-1,1),(-1,1)],[(-1,1),(-1,1)])
        [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
        sage: lattice_elements_in_box([[2,0],[0,2]],[(-2,2),(-1,1)],[(-1,1),(-1,1)])
        [(-1, 0), (0, 0), (1, 0)]
        sage: from hilbert_modgroup.all import *
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
        if norm_bound:
            norm = 1
        for i in range(n):
            alpha_i = 0.0
            for j in range(n):
                alpha_i = alpha_i + lattice_basis[i][j] * coordinate_vector[j]
            if alpha_i < lattice_bounds[i][0] or alpha_i > lattice_bounds[i][1]:
                # We need to discard this
                is_within_bounds = False
                break
            if norm_bound:
                norm = norm*alpha_i
        if norm_bound and abs(norm) > norm_bound:
            is_within_bounds = False
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
        sage: from hilbert_modgroup.all import *
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


cpdef find_candidate_cusps(p, z, use_lll=True, use_norm_bound=True, return_sigma_candidates=False,
                           initial_bd_d=None,
                           use_initial_bd_d=True):
    r"""
    Return candidates for closest cusp to the point ``z`` in the upper half-plane. 
    
    INPUT:
    
    - ``p`` -- object of type HilbertPullback
    - ``z`` -- point in the upper half-plane
    - ``use_lll`` -- boolean (default: `True`)  Use the LLL method to find a preliminary bounds
    - ``use_norm_bound`` -- boolean (default: `True`) Use the norm bound together ith the embedding bounds
    - ``return_sigma_candidates`` -- boolean (default: `False`) Return a list of sigma candidates only
    - ``initial_bd_d`` -- positive number (default: `None`) - an initial bound for the distance to nearest cusp.
    - ``use_initial_bd_d`` -- boolean (default: `False`) Use the initial bound
    
    OUTPUT:
    - list of tuples (sigma,rho) representing cusps 
    
    EXAMPLES::
    
        sage: from hilbert_modgroup.all import *
        sage: from hilbert_modgroup.pullback_cython import find_candidate_cusps
        sage: H1 = HilbertModularGroup(5)
        sage: P1 = HilbertPullback(H1)
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: find_candidate_cusps(P1,z)
        [(1, 0), (-1, -1), (0, 1), (1, -1)]
        sage: z=UpperHalfPlaneProductElement([CC(0,0.5),CC(0,1)])
        sage: find_candidate_cusps(P1,z)
        [(0, 1)]        
        sage: H2 = HilbertModularGroup(10)
        sage: P2 = HilbertPullback(H2)
        sage: z=UpperHalfPlaneProductElement([CC(0,1),CC(0,1)])
        sage: find_candidate_cusps(P2,z)
         [(1, 0),
         ...
         sage: len(_)
         10
        sage: z=UpperHalfPlaneProductElement([CC(2.58,0.5),CC(0.5,0.5)]) # long time (1 second)
        sage: len(find_candidate_cusps(P2,z)) # long time (1 second)
        271
        
    """
    ideal_basis = p._construct_ideal(1).integral_basis()
    lattice_basis = p.basis_matrix_ideal()
    n = len(lattice_basis[0])
    # Make lattice basis to a nested list to avoid creation of FreeModule elements
    lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
    ## Initial bound:
    if use_lll:
        candidate_cusp = p.get_heuristic_closest_cusp(z)
    else:
        candidate_cusp = None
    if candidate_cusp and use_initial_bd_d:
        dist = distance_to_cusp(p, candidate_cusp[0], candidate_cusp[1], z)
    else:
        dist = None
    if use_norm_bound:
        norm_bound_sigma = p._bound_for_sigma_norm(z, dist)
        norm_bound_sigma = norm_bound_sigma*(1+2**(-z[0].prec()/2)) # correct for numerical errors
    else:
        norm_bound_sigma = None
    if use_initial_bd_d and initial_bd_d and (not dist or dist > initial_bd_d):
        dist = initial_bd_d
        # norm_bound_rho = p._bound_for_rho_norm(z, dist)
    one = p.number_field()(1)
    coordinate_bounds = p._bound_for_sigma_coordinates(z, initial_bd_d=dist,use_initial_bd_d=use_initial_bd_d)
    embedding_bounds = p._bound_for_sigma_embeddings(z, initial_bd_d=dist,use_initial_bd_d=use_initial_bd_d)
    coordinate_bounds = [(-b,b) for b in coordinate_bounds]
    embedding_bounds = [(-b,b) for b in embedding_bounds]
    sigma_candidates_coordinates = lattice_elements_in_box(lattice_basis,
                                               embedding_bounds,
                                               coordinate_bounds, norm_bound=norm_bound_sigma)
    sigma_candidates = coordinates_to_ideal_elements(sigma_candidates_coordinates,
                                                             ideal_basis)
    if return_sigma_candidates:
        return sigma_candidates
    # To remove duplicates we keep track of the quotients rho/sigma (for sigma !=0)
    quotients = {}
    if not candidate_cusp or candidate_cusp[1] == 0:
        result = [(p.number_field()(1), p.number_field()(0))]
    elif candidate_cusp[0] == 0:
        result = [(p.number_field()(0), p.number_field()(1))]
    else:
        # We always have infinity since sigma=0 is always within the bounds.
        result = [candidate_cusp,(p.number_field()(1), p.number_field()(0))]
        q = candidate_cusp[0]/candidate_cusp[1]
        ngcd = gcd(candidate_cusp[0].norm(),candidate_cusp[1].norm())
        quotients = {q: (ngcd,(candidate_cusp[0],candidate_cusp[1]))}
    for s in sigma_candidates:
        if s == 0:
            continue
        rho_coordinate_bounds = p._bound_for_rho_coordinates(z, s, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d)
        rho_coordinate_bounds = [(-b,b) for b in rho_coordinate_bounds]
        rho_embedding_bounds = p._bound_for_rho_embeddings(z, s, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d)

        rho_candidates_coordinates = lattice_elements_in_box(lattice_basis,
                                                 rho_embedding_bounds,
                                                 rho_coordinate_bounds)
        rho_candidates = coordinates_to_ideal_elements(rho_candidates_coordinates,
                                                         ideal_basis)

        for r in rho_candidates:
            if r == 0 and (r,one) not in result:
                result.append((r,one))
            if r == 0:
                continue
            quo = r/s
            ngcd = gcd(r.norm(), s.norm())
            if quo in quotients:
                # We try to pick the representative cusp that has the smallest norm gcd
                if ngcd < quotients[quo][0]:
                    result.remove(quotients[quo][1])
                else:
                    continue
            result.append((r,s))
            quotients[quo] = (ngcd,(r,s))
    return result

cpdef find_closest_cusp(p, z, return_multiple=False, use_lll=True, use_norm_bound=True):
    r"""
    Return the closest cusp(s) to z
    
    INPUT:
    - ``p`` -- A HilbertPullback object
    - ``z`` -- point in the upper half-plane
    - ``return_multiple`` -- boolean (default False), set to True to return multiple closest cusp if it is no unique. 
    - ``use_lll`` -- boolean (default True) if True uses the LLL method to obtain an initial bound
    - ``use_norm_bound`` -- boolean (default True) if True uses the norm bound as well as the embedding bounds

    OUTPUT:
    - tuple representing a cusp (if return_multiple == False)
    - list of tuples representing cusps (if return_multiple == True) 
    
    EXAMPLES::

        sage: from hilbert_modgroup.all import *
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
        [(-a + 2, -2)]
        
    """
    cusp_candidates = find_candidate_cusps(p, z, use_lll=use_lll, use_norm_bound=use_norm_bound)
    min_cusp = cusp_candidates[0]
    min_d = distance_to_cusp(p,min_cusp[0], min_cusp[1], z)
    if return_multiple:
        min_cusp = [min_cusp]
    min_cusp_bound = p._bound_for_closest_cusp()
    eps = 2.0**(4-53) # Numerical precision for when we consider cusps to be equally close.
    # There may be many identical cusps in this list.
    # In particular there might be many
    for c in cusp_candidates[1:]:
        d = distance_to_cusp(p, c[0],c[1], z)
        if abs(d-min_d)<eps and return_multiple:
            min_cusp.append(c)
        if d < min_d - 2*eps:
            if return_multiple:
                min_cusp = [c]
            else:
                min_cusp = c
            min_d = d
        if d < min_cusp_bound:
            break
    # We have already filtered away identical cusps but as a final stage
    # we also check if the minimal cusps are equal to any of the fixed representatives.
    # other than infinity which is already taken care of
    if p.group().ncusps() == 1:
        return min_cusp
    result = []
    for cusp in p.group().cusps()[1:]:
        c, d = cusp.numerator(),cusp.denominator()
        quo = c/d
        if return_multiple:
            for r,s in min_cusp:
                if s == 0:
                    result.append((r, s))
                elif r == 0 and (0,1) not in result:
                    result.append((0,1))
                elif r == 0:
                    continue
                elif r/s == quo and (c,d) not in result:
                    result.append((c,d))
                elif r/s != quo:
                    result.append((r,s))
        elif min_cusp[1] != 0 and min_cusp[0]/min_cusp[1] == quo:
            result = (c,d)
        else:
            result = min_cusp
    return result

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
    
        sage: from hilbert_modgroup.all import *
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
    ideal_rs = r.parent().fractional_ideal(r,s)
    rlist = r.complex_embeddings()
    slist = s.complex_embeddings()
    n = len(slist)
    d = 1
    for i in range(n):
        d = d* ((-z._x[i] * slist[i] + rlist[i]) ** 2 * z._y[i] ** -1 + slist[i] ** 2 * z._y[i])
    d = ideal_rs.norm()**-1*d.sqrt()
    return d