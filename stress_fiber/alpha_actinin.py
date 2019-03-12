# encoding: utf-8
"""
αaas: α-actinin as a service
CDW 2019
"""

import itertools
import numpy as np

from ._spring import Spring
from . import _units
from . import _diffuse


class ActininHead:
    """One of the two heads of an alpha-actinin"""

    def __init__(self, actinin, side):
        self.actinin = actinin
        self.side = side  # Which side of the α-actinin is this on, 0 or 1
        self.binding_site = None  # None or g-actin site
        self._update_x()

    @property
    def bound(self):
        """Are you currently bound? Extrapolate from binding site link"""
        return self.binding_site is not None

    @property
    def nearest_binding_site(self):
        """The nearest actin binding site from neighboring tracts"""
        # Get all candidate actins
        actins = [t.actin for t in self.actinin.tract.neighbors]
        actins = list(itertools.chain(*actins))  # flatten
        # Find the binding site nearest our location
        near = [act.nearest(self.x) for act in actins]
        nearest = near[np.argmin(np.subtract([a.x for a in near], self.x))]
        return nearest

    @property
    def x(self):
        """Location derived from α-actinin"""
        # If bound, you are located at the binding site
        if self.bound:
            return self.binding_site.x
        # Else from parent α-actinin loc a distribution of parent lengths
        else:
            return self._x

    def _update_x(self):
        """Update the α-actinin-vibration-based estimate of x

        NOTE This is recalculating the x even when the actinin-head is bound.
        This isn't necessary and we could stop it by linking our updates to a
        global timestep clock.
        """
        spring = self.actinin.spring
        if self.side == 0:
            self._x = self.actinin.x - 0.5 * spring.bop_dx()
        elif self.side == 1:
            self._x = self.actinin.x + spring.rest + 0.5 * spring.bop_dx()

    def step(self):
        """Take a timestep: bind, unbind, or stay current"""
        self._update_x()
        if not self.bound:
            self._bind_or_not()
        else:
            self._unbind_or_not()

    def _bind_or_not(self):
        """Maybe bind? Can't say for sure."""
        actin_bs = self.nearest_binding_site
        rate = self._r12(abs(actin_bs.x - self.x))
        prob = rate * _units.constants.timestep
        if prob > np.random.rand():
            self.binding_site = actin_bs

    def _unbind_or_not(self):
        """Maybe unbind? Can't say for sure."""
        rate = self._r21()
        prob = rate * _units.constants.timestep
        if prob > np.random.rand():
            self.binding_site = None

    def _r12(self, dist):
        """Binding rate per second, given the distance to binding site.

        We take the binding rate to be dependent on the energy required to move
        from the current position to the binding site, with pre-exponential
        terms as in [1]_.

        .. [1] http://dx.doi.org/10.1098/rspb.2013.0697
        """
        tau = 72
        k = self.actinin.spring.k
        kT = _units.constants.kT
        rate = tau * np.exp(-(k * dist ** 2) / (2 * kT))
        return rate

    def _r21(self):
        """See _r12. We have a free energy release of ~ 8kcal/mol when α-actinin
        binds to vinculin [1]_. We'll use that as an approximation for the
        energy released when α-actinin binds to actin. This is equivalent to
        ~56 pN*nm per binding event. An estimated :math:`\Delta G_{12}` of
        56 pN*nm is approximately consistent with other estimates, such as [2]_
        which provides ~44 pN*nm. This is ~2x the energy released by the
        hydrolysis of ATP and so seems a reasonable estimate. We'll approximate
        the energy landscape as having a ~10% :math:`G_a` barrier of 6 pN*nm
        meaning the total energy barrier between the bound and unbound states is
        62 pN*nm.

        How much of that barrier the head has to surmount is dependent on how
        stretched the α-actinin backbone is at any point in time, thus
        :math:`\DeltaG_{a2}=62pN*nm-0.5*k*x^2`.

        We don't include a pre-exponential correction. Or, rather, assume it to
        be one.

        .. [1] https://dx.doi.org/10.1016/j.bpj.2012.08.044
        .. [2] https://dx.doi.org/10.1074/jbc.273.16.9570
        """
        deltaG = 62  # pN*nm energy barrier between unbound and bound
        U = self.actinin.energy
        A = 1
        kT = _units.constants.kT
        rate = A * np.exp(-(deltaG - U) / kT)
        return rate


class AlphaActinin:
    """A 1D α-actinin head

    As this head is considered as a semi-physical object, we care about its
    dimensions and material properties. The dimensions of α-actinin are readily
    available from crystallization, and more recently cryo-EM, studies.

    =========  =======  ======  ==========
    Property   Value    Units   Source
    =========  =======  ======  ==========
    Length     24-36    nm      [1]_, [2]_
    Thickness  0.5-6.5  nm      [1]_, [2]_
    =========  =======  ======  ==========

    Stiffness values are far harder to find. The flexural stiffness and
    persistence length of α-actinin have been measured or calculated via
    simulation multiple times, but I am unable to locate an extensional
    stiffness. So, let's calculate the Young's modulus from persistence length
    measurements made on generic α-helices and translate that to the specific
    dimensions of α-actinin. We take the persistence length of an α-helix to be
    80nm [3]_, and therefore for our tetrameric spectrin repeat backbone to be
    240nm. This assumes linear scaling, but is supported in [3]_ for scaling
    from single helices to coiled-coils. More structural detail can be found,
    amongst other places, here [4]_.

    Persistence length, :math:`L_p`, is related to the
    bending stiffness, :math:`B_s`, and therefore to Young's modulus, :math:`E`,
    via :math:`L_p=\frac{B_s}{k_B T}=\frac{E I}{k_B T}` where :math:`I` is the
    second moment of area. For a rod :math:`I=\frac{\pi r^4}{4}`.
    So our :math:`E` becomes :math:`E=\frac{4 L_p k_B T}{\pi r^4}`. What we
    really want is the spring constant, :math:`k` which we get from :math:`E` by
    :math:`k=\frac{E \pi r^2}{L_0}` and thus

    .. math:: k = \frac{4 L_p k_B T}{r^2 L_0}

    Where with :math:`L_p=240nm`, :math:`k_B=0.0138pN*nm/K`, :math:`T=288K`,
    :math:`r=1.5nm`, and :math:`L_0=36nm` to give :math:`k=3.75pN/nm`

    While we'll use this extensional stiffness below, it is worth noting that
    bending may be the primary mechanism by which α-actinin deforms [5]_.

    =========  =======  ======  ==========
    Property   Value    Units   Source
    =========  =======  ======  ==========
    Stiffness  3.75     pN/nm   [3]_
    =========  =======  ======  ==========

    .. [1] https://dx.doi.org/10.1007/s00018-008-8080-8
    .. [2] https://dx.doi.org/10.1016%2Fj.cell.2014.10.056
    .. [3] https://dx.doi.org/10.1103/PhysRevLett.97.248101
    .. [4] https://dx.doi.org/10.1371/journal.pone.0063633
    .. [5] https://doi.org/10.1016/S0092-8674(00)81980-7
    """

    def __init__(self, x, t):
        """An α-actinin head at location x in tract t"""
        self.heads = [ActininHead(self, 0), ActininHead(self, 1)]
        self.x = x
        self.tract = t
        self.spring = Spring(3.75, 36)  # See class doc for sources

    @property
    def bound(self):
        """Are you bound?"""
        return self.heads[0].bound or self.heads[1].bound

    @property
    def energy(self):
        """Current energy born by the stretched (or not) α-actinin"""
        if not (self.heads[0].bound and self.heads[1].bound):
            return 0
        dist = np.abs(self.heads[0].x - self.heads[1].x)
        return self.spring.energy(dist)

    def step(self):
        """Take a timestep"""
        if not self.bound:
            self.diffuse()
        [head.step() for head in self.heads]
        return

    def diffuse(self):
        """ Diffuse to a new location.
        We know the approximate dimensions of the α-actinin backbone are 24-36
        nm by 0.5-6.5 nm as specified in [Sjöblom_2008]_ and [Ribeiro_2014]_. We
        treat this as an ellipsoid of radii 18 by 3 nm to account for the heads
        at either end. 

        An alternate treatment would be to use the Stoke's radius for actinin
        described in [BNID:104395]

        .. [Sjöblom_2008] https://dx.doi.org/10.1007/s00018-008-8080-8
        .. [Ribeiro_2014] https://dx.doi.org/10.1016%2Fj.cell.2014.10.056
        .. [BNID:104395] https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104395
        """
        # NB This is checked well at this time, CDW20190221
        b, a = 18, 3
        f_drag = _diffuse.Drag.Ellipsoid.long_axis_translation(b, a)
        d_x = _diffuse.Dx(f_drag)
        self.x += d_x
        return d_x
