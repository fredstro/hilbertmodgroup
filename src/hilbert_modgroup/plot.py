r"""
Helper routine to make plots for the examples.


"""
from typing import ParamSpec

from sage.calculus.var import var
from sage.functions.other import ceil, floor
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.modules.free_module_element import vector
from sage.plot.all import Graphics, parametric_plot, point2d
from sage.plot.plot3d.implicit_plot3d import implicit_plot3d
from sage.plot.plot3d.shapes2 import point3d
from sage.repl.rich_output import get_display_manager
from sage.repl.rich_output.output_graphics3d import OutputSceneJmol as jmol
from sage.repl.rich_output.output_graphics3d import OutputSceneThreejs as threejs
from sage.rings.real_mpfr import RR

P = ParamSpec('P')


def plot_polygon(p: Polyhedron, basis: tuple, action: str = 'show', filename: str = '',
                 inside_polygon: bool = None,
                 norm_bound: float = 1.0,
                 curve_map: tuple = None,
                 norm_plot_factor: float = 1.0, **options: P.kwargs):
    """
    Plot a polygon in 2D or 3D.

    INPUT:

    - ``p`` -- polygon
    - ``basis`` -- basis for a lattice (usually given by an ideal)
    - ``action`` -- string (default 'show') can be 'show', 'save', 'return'
    - ``filename`` -- filename to save figure to if action = 'save'
    - `inside_bounds`` -- a list of bounds within which to restrict the lattice points
                            (to e.g. only plot inside the polygon)
    - ``norm_bound`` -- float (default: `None`) if set add a plot of xy= norm_bound
    - ``curve_map`` -- matrix - if set and norm_bound is set transform the curve xy=c with matrix
    - ``norm_plot_factor`` -- float (default: '1.0') - multiply the parameters of the xy=c plot
                                                       with this factor

    EXAMPLES::

        sage: from hilbert_modgroup.all import *
        sage: H1 = HilbertModularGroup(5)
        sage: K1 = H1.base_ring().number_field()
        sage: P1 = HilbertPullback(H1)
        sage: z = UpperHalfPlaneProductElement([0.1+I,0.2+0.5*I])
        sage: p=P1._candidate_integers_sigma(z,domain='preimage',return_polyhedron=True)
        sage: norm_args = {'norm_bound': 1.0}
        sage: plot_polygon(p,K1.fractional_ideal(1).integral_basis(),action='return',**norm_args)
        Graphics object consisting of 18 graphics primitives
        sage: x = ZZ['x'].gen()
        sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
        sage: H3=HilbertModularGroup(K3)
        sage: P3=HilbertPullback(H3)
        sage: z=UpperHalfPlaneProductElement([CC(0,1.0),CC(0,1.0),CC(0,1.0)])
        sage: p=P3._candidate_integers_sigma(z, domain='polytope', return_polyhedron=True)
        sage: norm_args = {'norm_bound': P3._bound_for_sigma_norm(z), 'viewer': 'jmol',
        ....:                   'curve_map': P3.basis_matrix_ideal(), 'norm_plot_factor': 1}
        sage: G = plot_polygon(p, [1,1,1], action='return', ticks=[2,2,2], polygon=False,
        ....: inside_polygon=True, xmin=-2, xmax=2, **norm_args)
        sage: type(G)
        <class 'sage.plot.plot3d.base.Graphics3dGroup'>
    """
    xmin = options.get('xmin') or floor(min(v[0] for v in p.vertices()))
    xmax = options.get('xmax') or ceil(max(v[0] for v in p.vertices()))
    ymin = options.get('ymin') or floor(min(v[1] for v in p.vertices()))
    ymax = options.get('ymax') or ceil(max(v[1] for v in p.vertices()))
    zmin = 0
    zmax = 1
    n = len(basis)
    if n == 3:
        zmin = options.get('zmin') or floor(min(v[2] for v in p.vertices()))
        zmax = options.get('zmax') or ceil(max(v[2] for v in p.vertices()))
    elif n > 3:
        raise NotImplementedError("This function is only implemented for degrees 2 and 3.")
    plot_options = {'point': False, 'line': 'black', 'fill': (0.9, 0.9, 0.9), 'wireframe': 'black',
                    'polygon': [0.9, 0.9, 0.9]}
    show_options = {'ticks': options.get('ticks', [xmax, ymax])}
    if n > 2:
        show_options['camera_position'] = options.get('camera_position', (2.3, 2.4, 2.0))
    if 'viewer' in options:
        show_options['viewer'] = options['viewer']
    # Remove unwanted options for the polot
    options = {k: v for k, v in options.items() if k in plot_options}
    if hasattr(basis[0], 'complex_embeddings'):
        labels = [r'$\varphi_1$', r'$\varphi_2$', r'$\varphi_3$'][0: n]
    else:
        labels = ['$X_1$', '$X_2$', '$X_3$'][0: n]

    plot_options.update(options)

    G = Graphics()
    G += p.plot(**plot_options)
    if n == 2:
        G = add_pts_degree2(G, basis, xmin, xmax, ymin, ymax, inside_polygon, p)
        G.set_axes_range(xmin - 1, xmax + 1, ymin - 1, ymax + 1)
        G.axes_labels(labels)

    else:
        G = add_pts_degree3(G, basis, xmin, xmax, ymin, ymax, zmin, zmax, inside_polygon,
                            curve_map, p, norm_bound)
    # Add also the norm bounds for degree 2
    if norm_bound and n == 2:
        B = curve_map or [[1, 0], [0, 1]]
        x_min = norm_bound / ymax / norm_plot_factor
        x_max = norm_plot_factor * xmax
        x = var('x')
        for (eps1, eps2) in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
            G += parametric_plot(((B[0][0] * eps1 * x + B[0][1] * eps2 * norm_bound / x),
                                  (B[1][0] * eps1 * x + B[1][1] * eps2 * norm_bound / x)),
                                 (x, x_min, x_max), zorder=1, color='black')
    if norm_bound and n == 3:
        x = var('x')
        y = var('y')
        z = var('z')
        B = curve_map or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        x1 = B[0][0] * x + B[0][1] * y + B[0][2] * z
        y1 = B[1][0] * x + B[1][1] * y + B[1][2] * z
        z1 = B[2][0] * x + B[2][1] * y + B[2][2] * z
        xmax = xmax*norm_plot_factor
        ymax = ymax*norm_plot_factor
        zmax = zmax*norm_plot_factor
        xmin = xmin * norm_plot_factor
        ymin = ymin * norm_plot_factor
        zmin = zmin * norm_plot_factor

        G += implicit_plot3d(abs(x1*y1*z1) == norm_bound, (x, xmin, xmax),
                             (y, ymin, ymax), (z, zmin, zmax), color='gray',
                             opacity=0.4, plot_points=100)
    return return_show(G, action, filename, **show_options)


def return_show(G: Graphics, action: str, filename: str, **show_options: P.kwargs):
    r"""
    Show, save, or return a graphics object based on the specified action.

    This function handles the output of a plot by either displaying it,
    saving it to a file, or returning it as a Graphics object.

    INPUT:

    - ``G`` -- a :class:`sage.plot.graphics.Graphics` object or a 3D graphics object
    - ``action`` -- string (default: ``'show'``); one of:

      - ``'show'`` -- display the plot using the appropriate backend
      - ``'save'`` -- save the plot to a file
      - ``'return'`` -- return the Graphics object without displaying

    - ``filename`` -- string; the filename to save the plot to (required if ``action='save'``)
    - ``**show_options`` -- additional keyword arguments passed to the ``show()`` or
      ``save()`` methods (e.g., ``ticks``, ``camera_position``, ``viewer``)

    OUTPUT:

    - If ``action='return'``, returns the Graphics object ``G``
    - If ``action='show'`` and no 3D viewer is available, returns the result of
      :func:`draw_3dplot_as_png`
    - Otherwise, returns ``None``

    EXAMPLES::

        sage: from sage.plot.all import Graphics, plot
        sage: from hilbert_modgroup.plot import return_show
        sage: G = plot(sin(x), (x, 0, 2*pi))
        sage: return_show(G, 'return', '') is G
        True

    ::

        sage: from sage.misc.temporary_file import tmp_filename
        sage: fn = tmp_filename(ext='.png')
        sage: return_show(G, 'save', fn)
        sage: import os; os.path.exists(fn)
        True

    Invalid action raises an error::

        sage: return_show(G, 'invalid', '')
        Traceback (most recent call last):
        ...
        ValueError: action: `invalid` is not a valid value.

    """
    dm = get_display_manager()
    show_3d = threejs in dm.supported_output() or jmol in dm.supported_output()
    if action == 'show' and hasattr(G, 'matplotlib'):
        G.matplotlib(**show_options)
    if action == 'show' and show_3d:
        G.show(**show_options)
    elif action == 'show':
        return draw_3dplot_as_png(G, **show_options)
    elif action == 'save' and filename:
        G.save(filename, **show_options)
    elif action == 'return':
        return G
    else:
        raise ValueError(f"action: `{action}` is not a valid value.")


def draw_3dplot_as_png(G: Graphics, **show_options: P.kwargs):
    r"""
    Render a 3D graphics object as a PNG image using matplotlib.

    This function is used as a fallback when no appropriate 3D viewer backend
    (such as threejs or jmol) is available. It saves the 3D plot to a temporary
    PNG file and displays it using matplotlib's ``imshow``.

    INPUT:

    - ``G`` -- a 3D graphics object (e.g., :class:`sage.plot.plot3d.base.Graphics3dGroup`)
    - ``**show_options`` -- keyword arguments passed to the ``save()`` method
      (e.g., ``ticks``, ``camera_position``)

    OUTPUT:

    This function does not return a value. It displays the plot using
    ``matplotlib.pyplot.imshow``.

    EXAMPLES:

    Create a simple 3D plot and render it as PNG::

        sage: from sage.plot.plot3d.shapes2 import sphere
        sage: from hilbert_modgroup.plot import draw_3dplot_as_png
        sage: G = sphere(center=(0,0,0), size=1, color='blue')
        sage: draw_3dplot_as_png(G, ticks=[1,1,1])  # renders to matplotlib

    Using with a Hilbert modular group polytope::

        sage: from hilbert_modgroup.all import *
        sage: x = ZZ['x'].gen()
        sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
        sage: H3 = HilbertModularGroup(K3)
        sage: P3 = HilbertPullback(H3)
        sage: z = UpperHalfPlaneProductElement([CC(0,1.0), CC(0,1.0), CC(0,1.0)])
        sage: p = P3._candidate_integers_sigma(z, domain='polytope', return_polyhedron=True)
        sage: norm_args = {'norm_bound': P3._bound_for_sigma_norm(z),
        ....:              'curve_map': P3.basis_matrix_ideal(), 'norm_plot_factor': 1}
        sage: G = plot_polygon(p, [1,1,1], action='return', ticks=[2,2,2], polygon=False,
        ....:                  inside_polygon=True, xmin=-2, xmax=2, **norm_args)
        sage: type(G)
        <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        sage: draw_3dplot_as_png(G, ticks=[2,2,2])  # renders the polytope as PNG
    """
    import matplotlib as mpl
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from sage.misc.temporary_file import tmp_filename
    mpl.rcParams['image.interpolation'] = 'bilinear'
    mpl.rcParams['image.resample'] = False
    mpl.rcParams['figure.figsize'] = [8.0, 6.0]
    mpl.rcParams['figure.dpi'] = 80
    mpl.rcParams['savefig.dpi'] = 100
    fn = tmp_filename(ext=".png")
    G.save(fn, **show_options)
    img = mpimg.imread(fn)
    show_options.pop('ticks')
    plt.imshow(img)


def add_pts_degree2(G: Graphics, basis: tuple, xmin: float, xmax: float, ymin: float, ymax: float,
                    inside_polygon: bool, p: Polyhedron) -> Graphics:
    """

        Add points to the plot G for the case of degree 2.

        INPUT:
        - ``G`` -- Graphics object to add points to
        - ``basis`` -- basis for a lattice (usually given by an ideal)
        - ``xmin`` -- minimum value of the x-axis
        - ``xmax`` -- maximum value of the x-axis
        - ``ymin`` -- minimum value of the y-axis
        - ``ymax`` -- maximum value of the y-axis
        - ``inside_polygon`` -- boolean, if True, only add points inside the polygon p
        - ``p`` -- polygon

        This function iterates over all integer coordinates within the given bounds,
        calculates the corresponding point in the lattice, and adds it to the plot
        if it is within the polygon (if specified).

        EXAMPLES::

            sage: from hilbert_modgroup.all import *
            sage: H1 = HilbertModularGroup(5)
            sage: K1 = H1.base_ring().number_field()
            sage: P1 = HilbertPullback(H1)
            sage: z = UpperHalfPlaneProductElement([0.1+I,0.2+0.5*I])
            sage: p = P1._candidate_integers_sigma(z, domain='preimage', return_polyhedron=True)
            sage: norm_args = {'norm_bound': 1.0}
            sage: G = Graphics()
            sage: basis = K1.fractional_ideal(1).integral_basis()
            sage: xmin = 0
            sage: xmax = 1
            sage: ymin = 0
            sage: ymax = 1
            sage: inside_polygon = False
            sage: from hilbert_modgroup.plot import add_pts_degree2
            sage: add_pts_degree2(G, basis, xmin, xmax, ymin, ymax, inside_polygon, p)
            Graphics object consisting of 2 graphics primitives

        """

    a, b = basis
    for i in range(floor(xmin), ceil(xmax) + 1):
        for j in range(floor(ymin), ceil(ymax) + 1):
            x = a * i + b * j
            if hasattr(x, 'complex_embeddings'):
                coordinates = [tuple(x.complex_embeddings())]
            else:
                coordinates = [(i, j)]
            if coordinates[0][0] < xmin or coordinates[0][0] > xmax:
                continue
            if coordinates[0][1] < ymin or coordinates[0][1] > ymax:
                continue

            if inside_polygon and coordinates[0] not in p:
                continue
            G += point2d(coordinates, zorder=1, color='black', size=20)
    return G


def add_pts_degree3(G, basis: tuple, xmin: float, xmax: float, ymin: float, ymax: float,
                    zmin: float, zmax: float, inside_polygon: bool = False,
                    curve_map: tuple = None, p: Polyhedron = None,
                    norm_bound: float = 1.0) -> Graphics:
    """
    Add points to the plot G for the case of degree 3

    INPUT:

    - ``basis`` - basis of the ideal
    - ``xmin``: float - minimum value of the x-axis
    - ``xmax``: float - maximum value of the x-axis
    - ``ymin``: float - minimum value of the y-axis
    - ``ymax``: float - maximum value of the y-axis
    - ``zmin``: float - minimum value of the z-axis
    - ``zmax``: float - maximum value of the z-axis
    - ``inside_polygon`` - boolean (default=False) if True, only add points inside the polygon p
    - ``curve_map`` - basis of the curve map
    - ``p`` - polygon
    - ``norm_bound`` - norm bound

    EXAMPLES:

            sage: from hilbert_modgroup.all import *
            sage: x = ZZ['x'].gen()
            sage: K3.<a> = NumberField(x^3-x^2-2*x+1)
            sage: H3=HilbertModularGroup(K3)
            sage: P3=HilbertPullback(H3)
            sage: z=UpperHalfPlaneProductElement([CC(0,1.0),CC(0,1.0),CC(0,1.0)])
            sage: p = P3._candidate_integers_sigma(z, domain='preimage', return_polyhedron=True)
            sage: norm_args = {'norm_bound': 1.0}
            sage: G = Graphics()
            sage: basis = K3.fractional_ideal(1).integral_basis()
            sage: xmin = 0
            sage: xmax = 1
            sage: ymin = 0
            sage: ymax = 1
            sage: zmin = 0
            sage: zmax = 1
            sage: inside_polygon = False
            sage: curve_map = None
            sage: from hilbert_modgroup.plot import add_pts_degree3
            sage: G = add_pts_degree3(G, basis, xmin, xmax, ymin, ymax, zmin, zmax)
            sage: type(G)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
    """
    a, b, c = basis
    coordinate_as_tuple = not hasattr(a, 'complex_embeddings')

    points_3d = [tuple((a * i + b * j + c * k).complex_embeddings())
                 if not coordinate_as_tuple else (i, j, k)
                 for i in range(floor(xmin), ceil(xmax) + 1)
                 for j in range(floor(ymin), ceil(ymax) + 1)
                 for k in range(floor(zmin), ceil(zmax) + 1)
                 ]
    for coordinates in points_3d:
        if coordinates[0] < xmin or coordinates[0] > xmax:
            continue
        if coordinates[1] < ymin or coordinates[1] > ymax:
            continue
        if coordinates[2] < zmin or coordinates[2] > zmax:
            continue

        if inside_polygon and coordinates not in p:
            continue
        if not curve_map:
            norm = coordinates[0] * coordinates[1] * coordinates[2]
        else:
            B = curve_map
            emb = B * vector(RR, coordinates)
            norm = emb[0] * emb[1] * emb[2]
        if norm_bound and abs(norm) > norm_bound + 0.01:
            G += point3d([coordinates], zorder=1, color='gray', size=30, opacity=0.75)
        else:
            G += point3d([coordinates], zorder=3, color='black', size=30)
    return G
