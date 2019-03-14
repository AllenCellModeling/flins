# encoding: utf-8
"""
Stagger around
"""

import numpy as np

from . import _units


class Drag:
    """Drag coefficients drawn from [Howard_2001]_, Pg 107 and [Berg_1983]_.
    We assume that our coefficient of viscosity, :math:`\eta`, defaults to that
    of water (0.0114 poise at 288K), referenced from _units

    .. [Howard_2001] http://books.google.com/books?vid=ISBN9780878933334
    .. [Berg_1983] https://press.princeton.edu/titles/112.html
    """

    class Cylinder:
        """Drag coefficients for cylinders.
        The cylinder has a length :math:`L`, a radius :math:`r`, and is in a
        fluid of viscosity :math:`\eta`.  This assumes :math:`L>>r` and an
        unbounded space. 
        """

        def long_axis_translation(L, r, eta=_units.constants.eta):
            """Drag coefficient for cylinder moving along long axis"""
            drag = 2 * np.pi * eta * L / (np.log(L / (2 * r)) - 0.20)
            return drag

        def short_axis_translation(L, r, eta=_units.constants.eta):
            """Drag coefficient for cylinder moving along the radial axis"""
            drag = 4 * np.pi * eta * L / (np.log(L / (2 * r)) + 0.84)
            return drag

        def long_axis_rotation(L, r, eta=_units.constants.eta):
            """Drag coefficient for cylinder rotating about the long axis"""
            drag = (1 / 3) * np.pi * eta * L ** 3 / (np.log(L / (2 * r)) - 0.66)
            return drag

        def short_axis_rotation(L, r, eta=_units.constants.eta):
            """Drag coefficient for cylinder rotating about the radial axis"""
            drag = 4 * np.pi * eta * r ** 2 * L
            return drag

    class Ellipsoid:
        """Drag coefficients for an ellipsoid.
        The ellipsoid has a minor axis radius of :math:`a`, a major axis radius
        of :math:`b`, and immersion in a fluid of viscosity :math:`\eta`. This
        assumes :math:`b>>a` and an unbounded space.
        """

        def long_axis_translation(b, a, eta=_units.constants.eta):
            """Drag coefficient for an ellipsoid moving along long axis"""
            drag = 4 * np.pi * eta * b / (np.log(2 * b / a) - 0.5)
            return drag

        def short_axis_translation(b, a, eta=_units.constants.eta):
            """Drag coefficient for an ellipsoid moving along the minor axis"""
            drag = 8 * np.pi * eta * b / (np.log(2 * b / a) + 0.5)
            return drag

        def long_axis_rotation(b, a, eta=_units.constants.eta):
            """Drag coefficient for an ellipsoid rotating about the long axis"""
            drag = (8 / 3) * np.pi * eta * b ** 3 / (np.log(2 * b / a) - 0.5)
            return drag

        def short_axis_rotation(b, a, eta=_units.constants.eta):
            """Drag coefficient for an ellipsoid rotating about the minor axis"""
            drag = 16 / 3 * np.pi * eta * a ** 2 * b
            return drag

    class Sphere:
        """Drag coefficients for a sphere.
        The sphere has a radius of :math:`r` and is embedded in a fluid of
        viscosity :math:`\eta`.
        """

        def translation(r, eta=_units.constants.eta):
            """Drag coefficient for a sphere moving through a fluid"""
            drag = 6 * np.pi * eta * r
            return drag

        def rotation(r, eta=_units.constants.eta):
            """Drag coefficient for a sphere rotating in a fluid"""
            drag = 8 * np.pi * eta * r ** 3
            return drag


def Dx(f_drag):
    """How far do we move because of diffusion subject to drag?
    From the Einstein-Smoluchowski relation via [Berg_1983]_ we know that the
    diffusion coefficient for a particle subject to a viscous drag, :math:`f`
    (in :math:`g/s`) is :math:`D=kT/f`. Further, from the same source and
    [Howard_2001]_ we know the drag on cylinders, ellipsoids, and spheres and
    calculate them above for reference. 

    We add a correction factor of 1/3.2 to account for the difference in
    diffusion between water and eukaryotic cytoplasm [Swaminathan_1997].

    An alternate treatment would be to use the Stoke's radius for all molecules,
    approximating them as spheres. There is an argument that anything beyond the
    spherical approximation is false precision. 

    .. [Berg_1983] https://press.princeton.edu/titles/112.html
    .. [Swaminathan_1997] https://doi.org/10.1016/S0006-3495(97)78835-0
    """
    # Gather info
    kT = _units.constants.kT
    D_corr = _units.world.D_cyto_corr  # Cytoplasm correction factor
    t = _units.world.timestep
    # Calculate
    D = kT / f_drag
    D *= D_corr
    std_dev = np.sqrt(2 * D * t)
    d_x = np.random.normal(0, std_dev)
    return d_x
