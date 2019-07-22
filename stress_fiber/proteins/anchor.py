# encoding: utf-8
"""
Don't stray far now.

Anchor keeps a linked binding site strongly attached to a given location. It is
intended to be used to represent focal adhesions and their movements. 
"""

from .base import Protein
from ..support.spring import Spring
from ..support.binding_site import BindingSite


class Anchor(Protein):
    """A strongly anchored spring that stays put
    We need to be able to specify the locations of, e.g., the ends of an actin
    being treated as bound to a focal adhesion or simply a boundary of the
    simulation. To do this while keeping our force solvers happy, we use very
    stiff springs attached to actin binding sites. 
    """

    def __init__(self, x, anchor_to=None, tract=None, k=1000, rest=0):
        """Create our anchor

        Parameters
        ----------
        x : float
            Initial x location of the anchor
        anchor_to : component with BindingSite, optional
            Protein segment with a binding site to anchor down
        tract : stress_fiber.space.Tract, optional
            Tract in which this anchor exists
        k : float, optional
            Stiffness of the anchor's attachment to its location in pN/nm [1000]
        rest : float, optional
            Rest length of the spring tethering the anchor in nm [0]
        """
        # Store tract if given, else create placeholders
        super().__init__("anchor", tract)
        # Store location and links
        self.x = x
        self.bs = BindingSite(self)
        if anchor_to is not None:
            self.bs.bind(anchor_to.bs)
        self.spring = Spring(k, rest)

    def __str__(self):
        """String representation of an anchor"""
        return "Anchor at %02f" % self.x

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

    def step(self):
        """Sit there and remain bound"""
        pass
