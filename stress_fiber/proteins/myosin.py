# encoding: utf-8
"""
Power up. 

Myosin provides motors to power active movement. Our myosin is a opposing dimer
with motors on both ends. Practically this means that we have two springs of
variable length connected by a rigid region. It looks like this::

    Head 1 bs    Connecting   Head 2 bs
         |        backbone         |     
         v            |            |     
                      v            |     
         *--/\/\--|                v     
                  |-------|              
             ^            |--/\/\--*     
             |                          
             |                ^         
             |                |         
          Head 1 /         Head 2 /
         S1 domain        S1 domain

"""

import math as m
import numpy as np

from .base import Protein
from ..base import Base
from ..support import units
from ..support import spring
from ..support import diffuse
from ..support import kinetics
from ..support import binding_site


class MyosinHead(Base):
    """A myosin head with binding site, base/connection, and kinetics. 
    Currently large parts of this are verbatim from gh:cdw/multifil.
    """

    def __init__(self, myosin, side):
        self.myosin = myosin
        self.side = side
        # Set up spring and state
        self.spring_properties = (
            {"rest": 5, "k": 5 / units.constants.kT},
            {"rest": 5, "k": 5 / units.constants.kT},
            {"rest": 0, "k": 5 / units.constants.kT},
        )
        self.state = 0  # start unbound
        self.bs = binding_site.BindingSite(self)

    def __str__(self):
        """String representation of myosin head"""
        state = ("isn't", "is weakly", "is tightly")[self.state]
        return "Myo head, %s bound" % state

    @property
    def state(self):
        """Options are 0,1,2 aka free, weakly bound, strongly bound."""
        return self._state

    @state.setter
    def state(self, new_state):
        """Update state and spring properties"""
        self._state = new_state
        sprops = self.spring_properties[new_state]
        try:
            self.spring.k = sprops["k"]
            self.spring.rest = sprops["rest"]
        except AttributeError:
            self.spring = spring.Spring(sprops["k"], sprops["rest"])

    @property
    def other_head(self):
        return self.myosin.heads[self.side ^ 1]

    @property
    def x_base(self):
        """The location of the base of the head"""
        return self.myosin.x_segments[1+self.side]

    @property
    def x_tip(self):
        """Return the current location"""
        # If bound, you are located at the linked site
        if self.bs.bound:
            return self.bs.linked.x
        # Else give a deviation from the parent myosin location
        if self.side == 0:
            return self.x_base - self._spring_length
        elif self.side == 1:
            return self.x_base + self._spring_length

    def force(self, x=None):
        """Force this myosin head exerts at tip
        """
        other = self.other_head
        if not other.bs.bound:
            return 0.0
        if x is None:
            x = self.x_tip
        total_length = abs(other.x_tip - x)
        total_rest = self.spring.rest + self.myosin.spring.rest + other.spring.rest
        k1, k2, k3 = self.spring.k, self.myosin.spring.k, other.spring.k
        total_k = (k1 * k2 * k3) / (k1 * k2 + k1 * k3 + k2 * k3)

    @property
    def _spring_length(self):
        """Current spring length, either from bound state or sampled"""
        if self.state == 0:
            spring_length = np.clip(self.spring.bop_length(), 0, np.Inf)
        else:
            spring_length = abs(self.x_base - self.bs.linked.x)
        return spring_length


class Myosin(Protein):
    def __init__(self, x, tract=None):
        """A two-sided myosin hexamer in a tract"""
        # Store tract if given, else create placeholders
        super().__init__("myosin", tract)
        # Create the motor heads at either end of the myosin
        self.heads = [MyosinHead(self, side) for side in (0, 1)]
        # Create spring backbone
        self.spring = spring.Spring(100, 40)  # FIXME These aren't validated
        # Create the locations we'll use to track this myosin
        self.springs = (self.heads[0].spring, self.spring, self.heads[1].spring)
        self.x = x  # the tip of head 1

    def __str__(self):
        """String rep of myosin"""
        b = str(self.boundaries)
        heads = str([h.state for h in self.heads])
        t = str(self.tract.loc)
        return "Myo at %s with heads in states %s on tract %s" % (b, heads, t)

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x
        self._update_x_segments()

    def _update_x_segments(self):
        """Update the segments of x
        _x_segments is [head 1 tip, head 1 base, head 2 base, head 2 tip]
        """
        self._sample_spring_lengths()
        seg_tots = np.cumsum(self._spring_lengths)
        self.x_segments = [self.x + seg_tot for seg_tot in seg_tots]
        self.x_segments.insert(0, self.x)

    def _sample_spring_lengths(self):
        """Bump the springs around a bit
        Note that the lengths are deterministic when the system is fully bound.
        I'd like these to vibrate.
        """
        if not self.fully_bound:
            self._spring_lengths = [np.clip(spring.bop_length(), 0, np.Inf)
                                    for spring in self.springs]
        else:
            # Solve for all lengths given total length, L and each spring's 
            # rest length rn, and stiffness, kn
            L = abs(self.heads[0].bs.linked.x - self.heads[1].bs.linked.x)
            (k1, r1), (k2, r2), (k3, r3) = [(s.k, s.rest) for s in self.springs]
            l1 = (k1 * k2 * r1 + k1 * k3 * r1 + k2 * k3 * (L - r2 - r3)) / (
                k1 * (k2 + k3) + k2 * k3
            )
            l2 = (k1 * k2 * r2 + k1 * k3 * (L - r1 - r3) + k2 * k3 * r2) / (
                k1 * (k2 + k3) + k2 * k3
            )
            l3 = (k1 * k2 * (L - r1 - r2) + k1 * k3 * r3 + k2 * k3 * r3) / (
                k1 * (k2 + k3) + k2 * k3
            )
            self._spring_lengths = [l1, l2, l3]

    @property
    def bound(self):
        """Bound at all?"""
        return any([h.bs.bound for h in self.heads])

    @property
    def fully_bound(self):
        """Bound fully?"""
        return all([h.bs.bound for h in self.heads])

    @property
    def boundaries(self):
        return self.x_segments[0], self.x_segments[3]


    def step(self):
        """Take a timestep"""
        self._update_x_segments()
        if not self.bound:
            self._freely_diffuse()
        [head.step() for head in self.heads]
        return


    def _freely_diffuse(self):
        """How myosin moves when unbound"""
        b, a = self.spring.rest, 15  # FIXME Myosin hexamer sizes not validated
        drag = diffuse.Drag.Ellipsoid.long_axis_translation(b, a)
        d_x = diffuse.Dx(drag)
        self.x += d_x
        if self.tract is not None:  # then stay in tract limits
            start, end = self.x_segments[0], self.x_segments[-1]
            limits = (0, self.tract.space.span)
            self.x, _ = diffuse.coerce_to_bounds(start, end, limits)
        return d_x
