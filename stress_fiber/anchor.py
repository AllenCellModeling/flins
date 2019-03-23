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
        self.binding_site = binding_site
        self.spring = Spring(k, rest)

    def force(self, x):
        """Force exerted by anchor if other end is stretched to x"""
        dist = abs(x - self.x)
        return self.Spring.force(dist)

    def energy(self, x):
        """Trick question, I don't want to consider this as storing energy.
        I suppose the anchor is dissipative, but really we don't want it
        contributing to calculations of the energy stored in cross linkers or
        motors. 
        """
        return 0
