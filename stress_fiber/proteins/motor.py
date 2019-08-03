# encoding: utf-8
r"""
Crank it. 

This is a simple simple motor. It isn't intended to do more than create some
sliding forces::

      Head 1        Head 2
         |             | 
         v             v
                           
         *--\/\/\/\/\--*
                  
                ^    
                |    
             Backbone

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


class Motor(Protein):
    def __init__(self, x, tract=None):
        # Store tract if given, else create placeholders
        super().__init__("motor", tract)
        self.x = x
        kT = units.constants.kT
        ks = (kT, kT, kT * 2)
        rs = (30.0, 30.0, 24.0)
        self.spring = [spring.Spring(k, r) for k, r in zip(ks, rs)]
        self.heads = [MotorHead(self, side) for side in (0, 1)]

    @property
    def state(self):
        return max([h.state for h in self.heads])

    @property
    def x(self):
        return self.locs[0]

    @x.setter
    def x(self, x):
        self._x = x

    def _update_x(self):
        """Ensure that x gets set each step depending on what's bound"""
        bound = self._which_bound
        if bound == "none":
            return  # x is tracked by motor
        elif bound == "left" or bound == "both":
            self.x = self.heads[0].bs.linked.x
        elif bound == "right":
            self.x = self.heads[1].bs.linked.x - self.spring[self.state].rest
        else:
            raise ValueError("_which_bound value not in expected set")

    @property
    def locs(self):
        """Location of each node in motor. 2 in this case"""
        bound = self._which_bound
        if bound == "none":
            x = self._x
            return (x, x + self.spring[self.state].rest)
        elif bound == "left":
            x = self.heads[0].bs.linked.x
            return (x, x + self.spring[self.state].rest)
        elif bound == "right":
            x = self.heads[1].bs.linked.x
            return (x - self.spring[self.state].rest, x)
        elif bound == "both":
            return [h.bs.linked.x for h in self.heads]
        else:
            raise ValueError("_which_bound value not in expected set")

    @property
    def bopped_locs(self):
        """Locations after perturbation, may not need this?"""
        # help us understand the bound-cases
        b1, b2 = [h.bs.bound for h in self.heads]
        none_bound = not b1 and not b2
        only_left_bound = b1 and not b2
        only_right_bound = not b1 and b2
        all_bound = b1 and b2
        # deal with the binding cases
        dx = self.spring[self.state].bop_dx()
        if none_bound:
            x = self.x
            return (x - 0.5 * dx, x + self.spring.rest + 0.5 * dx)
        elif only_left_bound:
            return (x, x + self.spring.rest + dx)
        elif only_right_bound:
            x = self.heads[1].bs.parent.x
            return (x - self.spring.rest - dx, x)
        elif all_bound:
            return [h.bs.parent.x for h in self.heads]

    @property
    def bound(self):
        return self.heads[0].bs.bound or self.heads[1].bs.bound

    @property
    def _which_bound(self):
        """Quick way to evaluate four cases: none, left, right, or both bound"""
        left, right = self.heads[0].bs.bound, self.heads[1].bs.bound
        if not left and not right:
            return "none"
        elif left and not right:
            return "left"
        elif not left and right:
            return "right"
        elif left and right:
            return "both"

    def _freely_diffuse(self):
        b, a = self.spring[0].rest, 15
        drag = diffuse.Drag.Ellipsoid.long_axis_translation(b, a)
        d_x = diffuse.Dx(drag)
        self.x += d_x
        if self.tract is not None:  # then stay in tract limits
            start, end = self.locs
            limits = (0, self.tract.space.span)
            self.x, _ = diffuse.coerce_to_bounds(start, end, limits)
        return d_x

    def step(self):
        if not self.bound:
            self._freely_diffuse()
        for head in self.heads:
            if head.state == 0:
                gact = self.tract.nearest_binding_site(head.x)
                if not gact.bs.bound:
                    head.step(bs=gact.bs)
            elif head.state == 1 or head.state == 2:
                length = np.abs(np.subtract(*self.locs))
                head.step(length=length)
            else:
                raise Exception("Head has state other than 0, 1, or 2")
        self._update_x()


class MotorHead(Base):
    """A head that tracks binding, rates, and states"""

    def __init__(self, motor, side):
        """ Create a motor head. 

        A note on states: typically we talk about a motor being in states that
        are one indexed (1, 2, ...) but since lists in python are zero indexed
        (0, 1, ...) we treat states as such to make it far easier to avoid
        off-by-one errors. 

        Parameters
        ==========
        motor: Motor
            The parent motor to this head
        side: 0 or 1
            Whether this head is on the left (0) or right (1) side of the motor
        """
        self.parent = motor
        self.side = side
        self.bs = binding_site.BindingSite(self)
        self.state = 0

    @property
    def other_head(self):
        return self.parent.heads[self.side ^ 1]

    @property
    def x(self):
        return self.parent.locs[self.side]

    def _r01(self, dist):
        """dist is distance to binding site"""
        return 100 * m.exp(-(0.25 * dist) ** 2)

    def _r10(self, length):
        """length is of motor backbone"""
        fe0 = self._free_energy(length, 0)
        fe1 = self._free_energy(length, 1)
        r01 = self._r01(
            abs(length - self.parent.spring[0].rest)
        )  # distance to bs is zero when bound
        return kinetics.reverse_rate(r01, fe0, fe1)

    def _r12(self, length):
        """See multifil.mh.Head for more documentation"""
        # fe2 = self._free_energy(dist, 1)
        # fe3 = self._free_energy(dist, 2)
        # scale, right_shift, speed = 0.6, 6.0, 0.2
        # return scale * (1 + m.tanh(right_shift + speed * (fe2 - fe3)))
        strain = self._strain(length, 1)
        scale, right_shift, speed = 48, 2.0, 0.5
        return scale * (1 - m.tanh(speed * (strain + right_shift))) + 4

    def _r21(self, length):
        fe1 = self._free_energy(length, 1)
        fe2 = self._free_energy(length, 2)
        r12 = self._r12(length)
        return kinetics.reverse_rate(r12, fe1, fe2)

    def _r20(self, length):
        strain = self._strain(length, 2)
        steepness, asymmetry, offset = 32, 3, 10
        return (m.sqrt(steepness * strain ** 2) - asymmetry * strain) + offset

    def _r02(self, dist):
        return 0.0

    def _strain(self, dist, state):
        return dist - self.parent.spring[state].rest

    def _free_energy(self, dist, state):
        """Free energy for given spring length and state"""
        kT = units.constants.kT
        if state == 0:
            return 0.0
        elif state == 1:
            return self.parent.spring[state].energy(dist) - 0.28 * 12 * kT
        elif state == 2:
            return self.parent.spring[state].energy(dist) - 0.68 * 12 * kT
        else:
            raise Exception("State other than 0, 1, or 2 given")

    def step(self, bs=None, length=None):
        """Take a timestep, transitioning to a new state as needed"""
        p = lambda r: kinetics.rate_to_prob(r, units.world.timestep)
        check = np.random.rand()
        if self.state == 0:
            if bs is None:
                raise Exception("Need binding site when in state 0")
            if self.other_head.bs.bound:
                other_fil = self.other_head.bs.linked.filament
                our_fil = bs.parent.filament
                if other_fil == our_fil:  # Don't self bind
                    return
            dist = abs(bs.parent.x - self.x)
            if p(self._r01(dist)) > check:
                self.bs.bind(bs)
                self.state = 1
        elif self.state == 1:
            if length is None:
                raise Exception("Need length when in state 1")
            if p(self._r12(length)) > check:
                self.state = 2
            elif 1 - p(self._r10(length)) < check:
                self.state = 0
                self.bs.unbind()
        elif self.state == 2:
            if length is None:
                raise Exception("Need length when in state 2")
            if p(self._r20(length)) > check:
                self.state = 0
                self.bs.unbind()
            elif 1 - p(self._r21(length)) < check:
                self.state = 1
        else:
            raise Exception("State other than 0, 1, or 2 given")

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
        """What force does this head exert or feel?"""
        force_fn = self.parent.spring[self.parent.state].force
        mult = -1 if self.x > self.other_head.x else 1
        return self._spring_property(force_fn, x) * mult

    def energy(self, x=None):
        """What energy is stored in the motor backbone?"""
        energy_fn = self.parent.spring[self.parent.state].energy
        return self._spring_property(energy_fn, x)
