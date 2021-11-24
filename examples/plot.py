r"""
Helper routine to make plots for the examples.


"""
from sage.all import floor, ceil, deepcopy, Graphics, point2d, point3d, var, parametric_plot,RR
from sage.modules.free_module_element import vector
from sage.plot.plot3d.implicit_plot3d import implicit_plot3d


def plot_polygon(p, basis, action='show', filename='', inside_polygon=None, norm_bound=None, curve_map=None,
                 norm_plot_factor=1.0, **options):
    """
    INPUT:
    - ``p`` -- polygon
    - ``basis`` -- basis for a lattice (usually given by an ideal)
    - ``action`` -- string (default 'show') can be 'show' or 'save'
    - ``filename`` -- filename to save figure to if action = 'save'
    - `inside_bounds`` -- a list of bounds within which to restrict the lattice points (to e.g. only plot inside the polygon)
    - ``norm_bound`` -- float (default: `None`) if set add a plot of xy= norm_bound
    - ``curve_map`` -- matrix - if set and norm_bound is set transform the curve xy=c with this matrix
    - ``norm_plot_factor`` -- float (default: '1.0') - extend the parameters of the xy=c plotwith this factor
    """
    xmin = options.get('xmin') or floor(min(v[0] for v in p.vertices()))
    xmax = options.get('xmax') or ceil(max(v[0] for v in p.vertices()))
    ymin = options.get('ymin') or floor(min(v[1] for v in p.vertices()))
    ymax = options.get('ymax') or ceil(max(v[1] for v in p.vertices()))
    n = len(basis)
    if n == 3:
        zmin = options.get('zmin') or floor(min(v[2] for v in p.vertices()))
        zmax = options.get('zmax') or ceil(max(v[2] for v in p.vertices()))
    elif n > 3:
        raise NotImplementedError("This plotting function is only implemented for degrees 2 and 3.")
    plot_options = {'point': False, 'line': 'black', 'fill': (0.9, 0.9, 0.9), 'wireframe': 'black',
                    'polygon': [0.9, 0.9, 0.9]}
    show_options = {'ticks': options.get('ticks', [xmax, ymax])}
    # Remove unwanted options for the polo
    for k in deepcopy(list(options)):
        if k not in plot_options:
            options.pop(k)
    if hasattr(basis[0], 'complex_embeddings'):
        labels = [r'$\varphi_1$', r'$\varphi_2$']
        if n == 3:
            labels += r'$\varphi_3$'
    else:
        labels = ['$X_1$', '$X_2$']
        if n == 3:
            labels += r'$X_3$'

    plot_options.update(options)

    G = Graphics();
    G += p.plot(**plot_options)
    if n == 2:
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
        G.set_axes_range(xmin - 1, xmax + 1, ymin - 1, ymax + 1)
        G.axes_labels(labels)

    else:
        a, b, c = basis
        for i in range(floor(xmin), ceil(xmax) + 1):
            for j in range(floor(ymin), ceil(ymax) + 1):
                for k in range(floor(zmin), ceil(zmax) + 1):
                    x = a * i + b * j + c * k
                    if hasattr(x, 'complex_embeddings'):
                        coordinates = [tuple(x.complex_embeddings())]
                    else:
                        coordinates = [(i, j, k)]
                    if coordinates[0][0] < xmin or coordinates[0][0] > xmax:
                        continue
                    if coordinates[0][1] < ymin or coordinates[0][1] > ymax:
                        continue
                    if coordinates[0][2] < zmin or coordinates[0][2] > zmax:
                        continue

                    if inside_polygon and coordinates[0] not in p:
                        continue
                    if not curve_map:
                        norm = coordinates[0][0]*coordinates[0][1]*coordinates[0][2]
                    else:
                        B = curve_map
                        emb = B*vector(RR,coordinates[0])
                        norm = emb[0]*emb[1]*emb[2]
                    if norm_bound and abs(norm) > norm_bound+0.01:
                        G += point3d(coordinates, zorder=1, color='gray', size=30,opacity=0.75)
                    else:
                        G += point3d(coordinates, zorder=3, color='black', size=30)

    # Add also the norm bounds for degree 2
    if norm_bound and n == 2:
        B = curve_map or [[1, 0], [0, 1]]
        x_min = norm_bound / ymax / norm_plot_factor
        x_max = norm_plot_factor * xmax
        x = var('x')
        for eps1 in [-1, 1]:
            for eps2 in [-1, 1]:
                G += parametric_plot(((B[0][0] * eps1 * x + B[0][1] * eps2 * norm_bound / x),
                                      (B[1][0] * eps1 * x + B[1][1] * eps2 * norm_bound / x)),
                                     (x, x_min, x_max), zorder=1, color='black')
    if norm_bound and n==3:
        x = var('x')
        y = var('y')
        z = var('z')
        B = curve_map or [[1,0,0],[0,1,0],[0,0,1]]
        x1 = B[0][0] * x + B[0][1] * y + B[0][2] * z
        y1 = B[1][0] * x + B[1][1] * y + B[1][2] * z
        z1 = B[2][0] * x + B[2][1] * y + B[2][2] * z
        xmax=xmax*norm_plot_factor
        ymax=ymax*norm_plot_factor
        zmax=zmax*norm_plot_factor
        xmin = xmin * norm_plot_factor
        ymin = ymin * norm_plot_factor
        zmin = zmin * norm_plot_factor

        G += implicit_plot3d(abs(x1*y1*z1) == norm_bound, (x, xmin, xmax), (y, ymin, ymax), (z, zmin, zmax), color='gray',
                             opacity=0.4,plot_points=100)
    if action == 'show':
        return G.show(**show_options)
    elif action == 'save' and filename:
        G.save(filename, **show_options)