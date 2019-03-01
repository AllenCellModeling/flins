# encoding: utf-8
"""
αaas: α-actinin as a service
CDW 2019
"""

from ._spring import Spring
from . import _units

import itertools
import numpy as np


class Actinin_head:
    """One of the two heads of an alpha-actinin"""

    def __init__(self, actinin, side):
        self.actinin = actinin
        self.side = side
        self.binding_site = None
        self.state = None  # None or g-actin site

    @property
    def x(self):
        """Location derived from actinin"""
        # If bound, you are located at the binding site
        if self.binding_site is not None:
            return self.binding_site.x
        # Else you have a default location from the α-actinin
        if side == 0:
            return self.actinin.x
        elif side == 1:
            return self.actinin.x + self.actinin.spring.rest

    def step(self, actins):
        """Take a timestep"""
        if self.state is None:
            self.bind_or_not()
        else:
            self.unbind_or_not()

    @property
    def nearest_binding_site(self):
        """The nearest actin binding site from neighboring tracts"""
        actins = [t.actin for t in self.actinin.tract.neighbors]
        actins = list(itertools.chain(actins))  # flatten
        near = [act.nearest(self.x) for act in actins]
        nearest = near[np.argmin(np.subtract([a.x for a in near], self.x))]
        return nearest

    def bind_or_not(self):
        """Maybe bind? Can't say for sure.
        TODO: This is very very wrong
        """
        actin_bs = self.nearest_binding_site
        rate = self._r12(actin_bs.x - self.x)
        if rate > np.random.rand():
            self.state = actin_bs
        return

    def unbind_or_not(self):
        """Maybe unbind? Can't say for sure.
        TODO: This is very very wrong
        """
        rate = self._r21(self.state.x - self.x)
        if rate > np.random.rand():
            self.state = None
        raise NotImplemented

    def _r12(self, dist):
        """ TODO Maybe split into a sep file with an eye to reuse
        TODO This is very very wrong
        """
        k = self.actinin.spring.k
        rest = self.actinin.spring.rest
        A = 2000
        ts = 1  # timestep
        rate = (
            float(A * np.sqrt(k / (2 * np.pi)) * np.exp(-.5 * k * (dist - rest) ** 2))
            * ts
        )
        return rate

    def _r21(self, dist):
        """See _r12"""
        # Predicate on a fraction of the energy of atp utilization
        # this is wrong and we need to extract it from lit values or take
        # educated guess and explain it
        g_atp = 13  # In units of RT
        atp = 5 * 10 ** -3
        adp = 30 * 10 ** -6
        phos = 3 * 10 ** -3
        deltaG = abs(-g_atp - log(atp / (adp * phos)))
        alphaDG = 0.28 * -deltaG
        etaDG = 0.68 * -deltaG
        r12 = self._r12(dist)
        try:
            rate = r12 / np.exp(alphaDG - etaDG)
        except ZeroDivisionError:
            rate = 1
        return float(rate)


class Alpha_actinin:
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

    =========  =======  ======  ==========
    Property   Value    Units   Source
    =========  =======  ======  ==========
    Stiffness  3.75     pN/nm   [3]_
    =========  =======  ======  ==========

    .. [1] https://dx.doi.org/10.1007/s00018-008-8080-8
    .. [2] https://dx.doi.org/10.1016%2Fj.cell.2014.10.056
    .. [3] https://dx.doi.org/10.1103/PhysRevLett.97.248101
    .. [4] https://dx.doi.org/10.1371/journal.pone.0063633
    """

    def __init__(self, x, t):
        """An α-actinin head at location x in tract t"""
        self.heads = [Actinin_head(self, 0), Actinin_head(self, 1)]
        self.x = x
        self.tract = t
        self.spring = Spring(3.75, 36)  # See class doc for sources

    @property
    def bound(self):
        """Are you bound?"""
        return any([h.state is not None for h in self.heads])

    def step(self):
        """Take a timestep"""
        if not self.bound:
            self.diffuse();
        [head.step() for head in self.heads]
        return

    def diffuse(self):
        """
        We know the approximate dimensions of the α-actinin backbone are 24-36
        nm by 0.5-6.5 nm as specified in [Sjöblom_2008]_ and [Ribeiro_2014]_. We
        treat this as an ellipsoid of radii 18 by 3 nm to account for the heads
        at either end. From the Einstein-Smoluchowski relation via [Berg_1983]_
        we have the diffusion coefficient for a particle subject to a viscous
        drag, :math:`f` as :math:`D=kT/f`. Further, from the same source we have
        the drag on an ellipsoidal particle along the long axis (:math:`a`) as
        :math:`f=\frac{4 \pi \eta a}{\ln 2a/b - 1/2}`. This gives us a diffusion
        coefficient of 

        .. math:: D = \frac{kT (\ln \frac{2a}{b} - \frac{1}{2})}{4 \pi \eta a}

        Where :math:`\eta` is the coefficient of viscosity for water, 0.0114
        poise at 288K. We add a correction factor of 1/3.2 to account for the
        difference in diffusion between water and eukaryotic cytoplasm
        [Swaminathan_1997].

        An alternate treatment would be to use the Stoke's radius for actinin
        described in [BNID:104395]

        .. [Sjöblom_2008] https://dx.doi.org/10.1007/s00018-008-8080-8
        .. [Ribeiro_2014] https://dx.doi.org/10.1016%2Fj.cell.2014.10.056
        .. [Berg_1983] https://press.princeton.edu/titles/112.html
        .. [Swaminathan_1997] https://doi.org/10.1016/S0006-3495(97)78835-0
        .. [BNID:104395] https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104395
        """
        # NOTE: This is checked well at this time, CDW20190221
        a, b = 18, 3
        eta = _units.constants.eta
        kT = _units.constants.boltzmann * _units.constants.temperature
        f_drag = (4*np.pi*eta*a)/(np.log(2*a/b) - 0.5)
        D = kT/f_drag  # Einstein-Smoluchowski relation for diffusion constant
        D *= 1/3.2  # Correction factor for diffusion in cytoplasm
        t = _units.constants.timestep
        std_dev = np.sqrt(2*D*t)
        d_x = np.random.normal(0, std_dev) 
        self.x += d_x
        return d_x
