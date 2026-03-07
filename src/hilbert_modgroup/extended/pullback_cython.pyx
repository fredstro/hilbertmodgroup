# %%cython
# cython: profile=True
"""
Cython versions of Extended Pullback algorithms.

Note: These are the algorithms that needs optimizing to make it all faster.
"""
from sage.all import ceil, floor, cartesian_product
from sage.arith.misc import gcd
from sage.structure.sage_object cimport SageObject
from hilbert_modgroup.upper_half_plane cimport UpperHalfPlaneProductElement__class
from sage.rings.number_field.number_field_element_base cimport NumberFieldElement_base

# Reuse shared functions from the parent package
from hilbert_modgroup.pullback_cython import (
    integral_coordinates_in_box,
    lattice_elements_in_box,
    coordinates_to_ideal_elements,
)

cpdef find_candidate_cusps(p, z, use_lll=True, use_norm_bound=True, return_sigma_candidates=False,
                           initial_bd_d=None,
                           use_initial_bd_d=True):
    r"""
    Return candidates for the closest cusp to the point ``z`` in the upper half-plane
    for an extended Hilbert modular group.

    This function enumerates cusp candidates `(\rho, \sigma)` by bounding
    the coordinates and embeddings of `\sigma` (and then `\rho` for each
    `\sigma`) using the geometry of the lattice ideal associated with ``p``.
    Duplicate cusps (those with the same quotient `\rho/\sigma`) are filtered,
    keeping the representative with the smallest norm gcd.

    INPUT:

    - ``p`` -- an :class:`~hilbert_modgroup.extended.pullback.ExtendedHilbertPullback` instance
    - ``z`` -- an element of :class:`~hilbert_modgroup.upper_half_plane.UpperHalfPlaneProductElement`
    - ``use_lll`` -- boolean (default: ``True``); if ``True``, use LLL reduction
      to find a preliminary heuristic closest cusp for tighter bounds
    - ``use_norm_bound`` -- boolean (default: ``True``); if ``True``, use the
      norm bound together with the embedding bounds when searching for sigma
      candidates
    - ``return_sigma_candidates`` -- boolean (default: ``False``); if ``True``,
      return only the list of sigma candidates without computing rho
    - ``initial_bd_d`` -- positive real number or ``None`` (default: ``None``);
      an initial upper bound for the distance to the nearest cusp
    - ``use_initial_bd_d`` -- boolean (default: ``True``); if ``True``, use the
      initial bound for the distance

    OUTPUT:

    A list of tuples `(\rho, \sigma)` representing candidate cusps
    `\rho/\sigma`.

    EXAMPLES::

        sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
        sage: from hilbert_modgroup.extended.pullback_cython import find_candidate_cusps
        sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
        sage: H1 = ExtendedHilbertModularGroup(5)
        sage: P1 = ExtendedHilbertPullback(H1)
        sage: z = UpperHalfPlaneProductElement([CC(0, 1), CC(0, 1)])
        sage: cusps = find_candidate_cusps(P1, z)
        sage: (P1.number_field()(1), P1.number_field()(0)) in cusps
        True
        sage: z = UpperHalfPlaneProductElement([CC(0, 0.5), CC(0, 1)])
        sage: find_candidate_cusps(P1, z)
        [(0, 1)]
        sage: sigmas = find_candidate_cusps(P1, z, return_sigma_candidates=True)
        sage: P1.number_field()(0) in sigmas
        True

    """
    ideal_basis = p._construct_ideal(1).integral_basis()
    lattice_basis = p.basis_matrix_ideal()
    n = len(lattice_basis[0])
    # Make lattice basis to a nested list to avoid creation of FreeModule elements
    lattice_basis = [[lattice_basis[i][j] for j in range(n)] for i in range(n)]
    # Initial bound:
    if use_lll:
        candidate_cusp = p.get_heuristic_closest_cusp(z)
    else:
        candidate_cusp = None
    if candidate_cusp and use_initial_bd_d:
        dist = distance_to_cusp_eg(p, candidate_cusp[0], candidate_cusp[1], z)
    else:
        dist = None
    if use_norm_bound:
        norm_bound_sigma = p._bound_for_sigma_norm(z, dist)
        norm_bound_sigma = norm_bound_sigma * (1 + 2 ** (-z[0].prec() / 2))  # correct for numerical errors
    else:
        norm_bound_sigma = None
    if use_initial_bd_d and initial_bd_d and (not dist or dist > initial_bd_d):
        dist = initial_bd_d
    one = p.number_field()(1)
    coordinate_bounds = p._bound_for_sigma_coordinates(z, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d)
    embedding_bounds = p._bound_for_sigma_embeddings(z, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d)
    coordinate_bounds = [(-b, b) for b in coordinate_bounds]
    embedding_bounds = [(-b, b) for b in embedding_bounds]
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
        result = [candidate_cusp, (p.number_field()(1), p.number_field()(0))]
        q = candidate_cusp[0] / candidate_cusp[1]
        ngcd = gcd(candidate_cusp[0].norm(), candidate_cusp[1].norm())
        quotients = {q: (ngcd, (candidate_cusp[0], candidate_cusp[1]))}
    for s in sigma_candidates:
        if s == 0:
            continue
        rho_coordinate_bounds = p._bound_for_rho_coordinates(z, s, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d)
        rho_coordinate_bounds = [(-b, b) for b in rho_coordinate_bounds]
        rho_embedding_bounds = p._bound_for_rho_embeddings(z, s, initial_bd_d=dist, use_initial_bd_d=use_initial_bd_d)

        rho_candidates_coordinates = lattice_elements_in_box(lattice_basis, rho_embedding_bounds, rho_coordinate_bounds)
        rho_candidates = coordinates_to_ideal_elements(rho_candidates_coordinates, ideal_basis)

        for r in rho_candidates:
            if r == 0 and (r, one) not in result:
                result.append((r, one))
            if r == 0:
                continue
            quo = r / s
            ngcd = gcd(r.norm(), s.norm())
            if quo in quotients:
                # We try to pick the representative cusp that has the smallest norm gcd
                if ngcd < quotients[quo][0]:
                    result.remove(quotients[quo][1])
                else:
                    continue
            result.append((r, s))
            quotients[quo] = (ngcd, (r, s))
    return result

cpdef find_closest_cusp(p, z, return_multiple=False, use_lll=True, use_norm_bound=True):
    r"""
    Return the closest cusp(s) to the point ``z`` in the upper half-plane
    for an extended Hilbert modular group.

    This function first obtains cusp candidates via
    :func:`find_candidate_cusps`, then evaluates the distance from ``z``
    to each candidate using :func:`distance_to_cusp_eg` and returns the
    one(s) with the smallest distance. When the group has multiple cusp
    classes, the result is normalized to use the fixed cusp representatives.

    INPUT:

    - ``p`` -- an :class:`~hilbert_modgroup.extended.pullback.ExtendedHilbertPullback` instance
    - ``z`` -- an element of :class:`~hilbert_modgroup.upper_half_plane.UpperHalfPlaneProductElement`
    - ``return_multiple`` -- boolean (default: ``False``); if ``True``, return
      a list of all cusps that are equally close (up to numerical precision)
    - ``use_lll`` -- boolean (default: ``True``); if ``True``, use LLL reduction
      for a preliminary heuristic bound
    - ``use_norm_bound`` -- boolean (default: ``True``); if ``True``, use the
      norm bound together with embedding bounds

    OUTPUT:

    - A tuple `(\rho, \sigma)` representing the closest cusp (if
      ``return_multiple`` is ``False``).
    - A list of tuples representing all equally-close cusps (if
      ``return_multiple`` is ``True``).

    EXAMPLES::

        sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
        sage: from hilbert_modgroup.extended.pullback_cython import find_closest_cusp
        sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
        sage: H1 = ExtendedHilbertModularGroup(5)
        sage: P1 = ExtendedHilbertPullback(H1)
        sage: z = UpperHalfPlaneProductElement([CC(0, 1), CC(0, 1)])
        sage: find_closest_cusp(P1, z)
        (1, 0)
        sage: find_closest_cusp(P1, z, return_multiple=True)
        [(1, 0), (0, 1)]
        sage: z = UpperHalfPlaneProductElement([CC(0, 0.5), CC(0, 1)])
        sage: find_closest_cusp(P1, z)
        (0, 1)
        sage: find_closest_cusp(P1, z, return_multiple=True)
        [(0, 1)]

    """
    cusp_candidates = find_candidate_cusps(p, z, use_lll=use_lll, use_norm_bound=use_norm_bound)
    min_cusp = cusp_candidates[0]
    min_d = distance_to_cusp_eg(p, min_cusp[0], min_cusp[1], z)
    if return_multiple:
        min_cusp = [min_cusp]
    min_cusp_bound = p._bound_for_closest_cusp()
    eps = 2.0 ** (4 - 53)  # Numerical precision for when we consider cusps to be equally close.
    # There may be many identical cusps in this list.
    for c in cusp_candidates[1:]:
        d = distance_to_cusp_eg(p, c[0], c[1], z)
        if abs(d - min_d) < eps and return_multiple:
            min_cusp.append(c)
        if d < min_d - 2 * eps:
            if return_multiple:
                min_cusp = [c]
            else:
                min_cusp = c
            min_d = d
        if d < min_cusp_bound:
            break
    # We have already filtered away identical cusps but as a final stage
    # we also check if the minimal cusps are equal to any of the fixed representatives.
    if p.group().ncusps() == 1:
        return min_cusp
    result = []
    for cusp in p.group().cusps()[1:]:
        c, d = cusp.numerator(), cusp.denominator()
        quo = c / d
        if return_multiple:
            for r, s in min_cusp:
                if s == 0:
                    result.append((r, s))
                elif r == 0 and (0, 1) not in result:
                    result.append((0, 1))
                elif r == 0:
                    continue
                elif r / s == quo and (c, d) not in result:
                    result.append((c, d))
                elif r / s != quo:
                    result.append((r, s))
        elif min_cusp[1] != 0 and min_cusp[0] / min_cusp[1] == quo:
            result = (c, d)
        else:
            result = min_cusp
    return result

cpdef distance_to_cusp_eg(SageObject p, NumberFieldElement_base r,
                          NumberFieldElement_base s,
                          UpperHalfPlaneProductElement__class z):
    r"""
    Return the distance from a point in the upper half-plane product to a cusp
    `\rho/\sigma` for an extended Hilbert modular group.

    The distance is computed as

    .. MATH::

        d(z, \rho/\sigma) = N(\mathfrak{a})^{-1}
        \prod_{i=1}^{n} \sqrt{(-x_i \sigma^{(i)} + \rho^{(i)})^2 y_i^{-1}
        + (\sigma^{(i)})^2 y_i}

    where `z = (x_1 + i y_1, \ldots, x_n + i y_n)`, the superscript
    `{}^{(i)}` denotes the `i`-th complex embedding, and `\mathfrak{a}` is
    the ideal generated by `(\rho, \sigma)` with respect to the group.

    INPUT:

    - ``p`` -- an :class:`~hilbert_modgroup.extended.pullback.ExtendedHilbertPullback`
      instance (typed as ``SageObject`` for Cython compatibility)
    - ``r`` -- a number field element; the numerator of the cusp
    - ``s`` -- a number field element; the denominator of the cusp
    - ``z`` -- an element of :class:`~hilbert_modgroup.upper_half_plane.UpperHalfPlaneProductElement`

    OUTPUT:

    A positive real number giving the distance from ``z`` to the cusp
    `\rho/\sigma`.

    EXAMPLES::

        sage: from hilbert_modgroup.extended.all import ExtendedHilbertModularGroup, ExtendedHilbertPullback
        sage: from hilbert_modgroup.extended.pullback_cython import distance_to_cusp_eg
        sage: from hilbert_modgroup.upper_half_plane import UpperHalfPlaneProductElement
        sage: H1 = ExtendedHilbertModularGroup(5)
        sage: P1 = ExtendedHilbertPullback(H1)
        sage: K = P1.number_field()
        sage: z = UpperHalfPlaneProductElement([CC(0, 1), CC(0, 1)])
        sage: distance_to_cusp_eg(P1, K(0), K(1), z)  # abs tol 1e-10
        1.0
        sage: distance_to_cusp_eg(P1, K(1), K(1), z)  # abs tol 1e-10
        2.0
        sage: z = UpperHalfPlaneProductElement([CC(0, 0.5), CC(0, 1)])
        sage: distance_to_cusp_eg(P1, K(0), K(1), z)  # abs tol 1e-10
        0.707106781186548

    """
    cdef list rlist, slist
    ideal_rs = p.group().ideal((r, s))
    rlist = r.complex_embeddings()
    slist = s.complex_embeddings()
    n = len(slist)
    d = 1
    for i in range(n):
        d = d * ((-z._x[i] * slist[i] + rlist[i]) ** 2 * z._y[i] ** -1 + slist[i] ** 2 * z._y[i])
    d = ideal_rs.norm() ** -1 * d.sqrt()
    return d
