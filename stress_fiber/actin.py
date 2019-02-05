# encoding: utf-8
"""
Actin' up
CDW 2019
"""

import numpy as np

class Actin:
    """An actin filament, in all its glory. 
    """
    def __init__(self, loc, act_pair):
        """Create an actin filament
        Parameters
        ----------
        loc: float
            The x location of our actin fil.
        act_pair: int
            The number of actin pairs making up the filament, to start. If
            positive then the filament is facing to the right, if negative then
            facing to the left. 
        """
        self.loc = loc
        self.act_pair = act_pair
        # Calculate actin rise and run, store for later use
        # Numbers derived from Howard, 2001, page 125
        monomers_per_polymer = 26
        polymer_base_length = 72.0
        polymer_base_turns = 12.0
        rev = 2*np.pi    
        self._pitch = polymer_base_turns * rev / monomers_per_polymer
        self._rise = polymer_base_length / monomers_per_polymer

    @property
    def length(self):
        """How long are you in nm?"""
        return abs(self.act_pair * self._rise)

    @property
    def bounds(self):
        """Where do you go from and to?"""
        return (self.loc, self.loc + self.act_pair * self._rise)

    def lengthen_shorten(self):
        """TODO, also to derive from Howard, 2001"""
        raise NotImplementedError

