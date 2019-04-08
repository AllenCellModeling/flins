# encoding: utf-8
"""
Don't stray far now
CDW 2019
"""

import numpy as np

from . import _spring
from . import _binding_site


class Anchor:
    """A strongly anchored spring that stays put
    We need to be able to specify the locations of, e.g., the ends of an actin
    being treated as bound to a focal adhesion or simply a boundary of the
    simulation. To do this while keeping our force solvers happy, we use very
    stiff springs attached to actin binding sites. 
    """

    def __init__(self, x, bind_to=None, k=100, rest=10):
        self.x = x
        self.bs = _binding_site.BindingSite(self)
        if binding_site is not None:
            self.bs.bind(bind_to)
        self.spring = _spring.Spring(k, rest)

    def __str__(self):
        """String representation of an anchor"""
        return "Anchor at %02f"%self.x
        
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
