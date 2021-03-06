# encoding: utf-8
"""
αaas: α-actinin as a service

But seriously, α-actinin is the primary cross-linker in this system. It requires
two heads separated by a flexible (or extensible) backbone that are able to bind
adjacent actin filaments.
"""

import numpy as np

from .base import Protein, Head
from ..support import spring
from ..support import units
from ..support import diffuse
from ..support import kinetics


class AlphaActinin(Protein):
    r"""A 1D α-actinin molecule with heads on either end

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
    amongst other places, in Autore_2013_.

    Persistence length, :math:`L_p`, is related to the bending stiffness,
    :math:`B_s`, and therefore to Young's modulus, :math:`E`, via
    :math:`L_p=\frac{B_s}{k_B T}=\frac{E I}{k_B T}` where :math:`I` is the
    second moment of area. For a rod :math:`I=\frac{\pi r^4}{4}`.  So our
    :math:`E` becomes :math:`E=\frac{4 L_p k_B T}{\pi r^4}`. What we really
    want is the spring constant, :math:`k` which we get from :math:`E` by
    :math:`k=\frac{E \pi r^2}{L_0}` and thus

    .. math:: k = \frac{4 L_p k_B T}{r^2 L_0}

    Where with :math:`L_p=240nm`, :math:`k_B=0.0138pN*nm/K`, :math:`T=288K`,
    :math:`r=1.5nm`, and :math:`L_0=36nm` to give :math:`k=3.75pN/nm`

    While we'll use this extensional stiffness below, it is worth noting that
    bending may be the primary mechanism by which α-actinin deforms (Grum_1999_).

    =========  =======  ======  ==========
    Property   Value    Units   Source
    =========  =======  ======  ==========
    Stiffness  3.75     pN/nm   [3]_
    =========  =======  ======  ==========

    .. [1] https://dx.doi.org/10.1007/s00018-008-8080-8
    .. [2] https://dx.doi.org/10.1016%2Fj.cell.2014.10.056
    .. [3] https://dx.doi.org/10.1103/PhysRevLett.97.248101
    .. _Autore_2013: https://dx.doi.org/10.1371/journal.pone.0063633
    .. _Grum_1999: https://doi.org/10.1016/S0092-8674(00)81980-7
    """

    def __init__(self, x, tract=None):
        """An α-actinin at location x in a tract"""
        # Store tract if given, else create placeholders
        super().__init__("actinin", tract)
        # Create the spring that is our actinin and remember passed values
        self.spring = spring.Spring(3.75, 36)  # See class doc for sources
        self.x = x
        # Create the heads at either end of the actinin
        self.heads = [ActininHead(self, 0), ActininHead(self, 1)]

    def __str__(self):
        """String representation of α-actinin"""
        tract_str = "NA" if self.tract is None else str(self.tract.loc)
        x_str = "%.1f" % self.x
        bound_str = "%i bound heads" % sum([h.bs.bound for h in self.heads])
        return "α-act with %s at x=%s on tract %s" % (bound_str, x_str, tract_str)

    @property
    def bound(self):
        """Are you bound?"""
        return self.heads[0].bs.bound or self.heads[1].bs.bound

    @property
    def fully_bound(self):
        """Are both sides bound?"""
        return self.heads[0].bs.bound and self.heads[1].bs.bound

    @property
    def energy(self):
        """Current energy born by the stretched (or not) α-actinin"""
        if not (self.heads[0].bs.bound and self.heads[1].bs.bound):
            return 0
        dist = np.abs(self.heads[0].x - self.heads[1].x)
        return self.spring.energy(dist)

    def step(self):
        """Take a timestep"""
        if not self.bound:
            self.freely_diffuse()
        else:
            if self.heads[0].bs.bound:
                self.x = self.heads[0].bs.linked.x
            else:
                self.x = self.heads[1].bs.linked.x - self.spring.rest
        [head.step() for head in self.heads]
        return

    def freely_diffuse(self):
        """ Diffuse to a new location.
        We know the approximate dimensions of the α-actinin backbone are 24-36
        nm by 0.5-6.5 nm as specified in Sjöblom_2008_ and Ribeiro_2014_. We
        treat this as an ellipsoid of radii 18 by 3 nm to account for the heads
        at either end.

        An alternate treatment would be to use the Stoke's radius for actinin
        described in BNID_104395_.

        .. _Sjöblom_2008: https://dx.doi.org/10.1007/s00018-008-8080-8
        .. _Ribeiro_2014: https://dx.doi.org/10.1016%2Fj.cell.2014.10.056
        .. _BNID_104395: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104395
        """
        # NB This is checked well at this time, CDW20190221
        b, a = 18, 3
        f_drag = diffuse.Drag.Ellipsoid.long_axis_translation(b, a)
        d_x = diffuse.Dx(f_drag)
        self.x += d_x
        if self.tract is not None:  # then derive diffusion limits from tract
            start, end = self.x, self.x + self.spring.rest
            space_limits = (0, self.tract.space.span)
            self.x, _ = diffuse.coerce_to_bounds(start, end, space_limits)
        return d_x


class ActininHead(Head):
    """One of the two heads of an α-actinin"""

    def __init__(self, actinin, side):
        super().__init__(actinin, side)
        self.address = (actinin.address[:], ("actininhead", side))
        self._update_x()

    def __str__(self):
        """String representation of α-actinin head"""
        x_str = "%0.1f" % self.x
        state_str = "unbound" if not self.bs.bound else "bound"
        return "α-act head %s at x=%s" % (state_str, x_str)

    @property
    def x(self):
        """Location derived from α-actinin"""
        # If bound, you are located at the linked site
        if self.bs.bound:
            return self.bs.linked.x
        # Else from parent α-actinin loc a distribution of parent lengths
        else:
            return self._x

    def _update_x(self):
        """Update the α-actinin-vibration-based estimate of x

        NOTE This is recalculating the x even when the actinin-head is bound.
        This isn't necessary and we could stop it by linking our updates to a
        global timestep clock.
        """
        spring = self.parent.spring
        if self.side == 0:
            self._x = self.parent.x - 0.5 * spring.bop_dx()
        elif self.side == 1:
            self._x = self.parent.x + spring.rest + 0.5 * spring.bop_dx()

    def step(self):
        """Take a timestep: bind, unbind, or stay current"""
        self._update_x()
        if not self.bs.bound:
            self._bind_or_not()
        else:
            self._unbind_or_not()

    def _bind_or_not(self):
        """Maybe bind? Can't say for sure."""
        gactin = self.parent.tract.nearest_binding_site(self.x)
        if gactin.bs.bound:  # don't bind if site is already taken
            return
        if self.other_head.bs.bound:
            if self.other_head.bs.linked.filament == gactin.filament:  # don't self bind
                return
        rate = self._r01(abs(gactin.x - self.x))
        prob = kinetics.rate_to_prob(rate, units.world.timestep)
        if prob > np.random.rand():
            self.bs.bind(gactin.bs)

    def _unbind_or_not(self):
        """Maybe unbind? Can't say for sure."""
        rate = self._r10()
        prob = kinetics.rate_to_prob(rate, units.world.timestep)
        if prob > np.random.rand():
            self.bs.unbind()

    def _r01(self, dist):
        """Binding rate per second, given the distance to binding site.

        We take the binding rate to be dependent on the energy required to move
        from the current position to the binding site, with pre-exponential
        terms as in [1]_.

        .. [1] http://dx.doi.org/10.1098/rspb.2013.0697
        """
        tau = 72
        k = self.parent.spring.k
        kT = units.constants.kT
        rate = tau * np.exp(-(k * dist ** 2) / (2 * kT))
        return rate

    def _r10(self):
        r"""See _r01. We have a free energy release of ~ 8kcal/mol when α-actinin
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
        U = self.parent.energy
        A = 1
        kT = units.constants.kT
        rate = A * np.exp(-(deltaG - U) / kT)
        return rate

    def _spring_property(self, spring_prop, x=None):
        """Calculate a spring property, energy or force"""
        # If other head isn't bound, can't bear energy/strain
        if not self.other_head.bs.bound:
            return 0
        # If no x is given, use current location
        if x is None:
            x = self.x
        length = abs(x - self.other_head.x)
        return spring_prop(length)

    def force(self, x=None):
        r"""What force does this α-actinin head exert or feel?

        This accounts for the fact that heads feel equal and opposite forces
        and that these forces change as the spring is compressed or extended
        from rest. Let's think of two heads, `a` and `b` at either end of an
        α-actinin spring. The default force sign returned from the spring is
        negative when the spring is shortened and positive when it is
        lengthened. The desired sign is that which reflects the force exerted by
        the spring on the binding site.

        ====  ===========  ======  =======  =======  ============
               System state              Force direction
        -------------------------  ------------------------------
        Head  Orientation  Spring  Default  Desired  Flip needed?
        ====  ===========  ======  =======  =======  ============
         A        A>B      Short     \-       \+         Yes
         A        A>B      Long      \+       \-         Yes
         A        B>A      Short     \-       \-         No
         A        B>A      Long      \+       \+         No
         B        A>B      Short     \-       \-         No
         B        A>B      Long      \+       \+         No
         B        B>A      Short     \-       \+         Yes
         B        B>A      Long      \+       \-         Yes
        ====  ===========  ======  =======  =======  ============

        Looking at this it becomes obvious that force direction flips are needed
        in cases where the current head is the right-most of the two.
        """
        force_fn = self.parent.spring.force
        mult = -1 if self.x > self.other_head.x else 1
        return self._spring_property(force_fn, x) * mult

    def energy(self, x=None):
        """What energy is stored in the α-actinin backbone?
        This exists because we want to be able to propose alternate ActininHead
        locations and find the energy without changing states.
        """
        energy_fn = self.parent.spring.energy
        return self._spring_property(energy_fn, x)
