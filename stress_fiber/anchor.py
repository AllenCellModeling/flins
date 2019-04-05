# encoding: utf-8
"""
Don't stray far now
CDW 2019
"""

import numpy as np

from ._spring import Spring


class Anchor:
    """A strongly anchored spring that stays put
    We need to be able to specify the locations of, e.g., the ends of an actin
    being treated as bound to a focal adhesion or simply a boundary of the
    simulation. To do this while keeping our force solvers happy, we use very
    stiff springs attached to actin binding sites. 
    """

    def __init__(self, x, binding_site=None, k=100, rest=10):
        self.x = x
        self._binding_site = binding_site
        self.spring = Spring(k, rest)

    def __str__(self):
        """String representation of an anchor"""
        return "Anchor at %02f"%self.x

    @property
    def binding_site(self):
        """Site that anchor is bound to, or none"""
        return self._binding_site

    @binding_site.setter
    def binding_site(self, bs):
        if self._binding_site is not None:
            self._binding_site.link = None
        if bs is not None:
            bs.link = self
        self._binding_site = bs
        
    def force(self, x):
        """Force exerted by anchor if other end is stretched to x"""
        dist = abs(x - self.x)
        return self.spring.force(dist)

    def energy(self, x):
        """Trick question, I don't want to consider this as storing energy.
        I suppose the anchor is dissipative, but really we don't want it
        contributing to calculations of the energy stored in cross linkers or
        motors. 
        """
        # TODO this is used to calculate actin boping-around and as such needs
        # to be passed back relatively high
        return 0
