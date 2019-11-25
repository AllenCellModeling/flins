# encoding: utf-8
"""
Just a slice of the cake^h^h^h^h space.

Tracts are sections of space that extend in a long direction and have adjacent
neighboring tracts that also extend in a long direction.
"""

import copy
import itertools
import numpy as np

from ..base import Base
from ..support import names


class Tract(Base):
    """A single tract in a Space"""

    def __init__(self, loc, space):
        """A single tract in a space"""
        self.loc = loc
        self.space = space
        self._neighbors = None
        self._reachable = None
        self.mols = {}
        self.mols_named = {}
        self.address = (("tract", loc),)

    def __str__(self):
        """String representation of tract"""
        loc_str = str(self.loc)
        mols_counts = ["%i %ss" % (len(v), k) for k, v in self.mols.items()]
        mols_str = ", ".join(mols_counts)
        return "Tract %s with %s" % (loc_str, mols_str)

    @property
    def neighbors(self):
        """Who are your neighbors?
        We don't want to calculate this each time, but we can't populate it on
        creation as there are neighbors that haven't yet been created. So we
        populate this list on the first call and then reference the stored
        version thereafter."""
        if self._neighbors is None:
            if self.space is None:
                self._neighbors = None
            else:
                self._neighbors = self.space.neighbors(self.loc)
        return self._neighbors

    @property
    def reachable(self):
        """Reachable tracts from here. Neighbors and self.

        See neighbors documentation for creation method.
        """
        if self._reachable is None:
            if self.space is None:
                self._reachable = [self]
            else:
                self._reachable = copy.copy(self.neighbors)
                self._reachable.append(self)
        return self._reachable

    def add_mol(self, kind, mol):
        """Add a molecule to our lists and dicts thereof"""
        # Create entries for this type of protein if seeing for first time
        if kind not in self.mols:
            self.mols[kind] = []
        if kind not in self.mols_named:
            self.mols_named[kind] = {}
        # Make sure we aren't already tracking this protein
        assert mol not in self.mols[kind], "Can't add mol to tract twice"
        # Create id for this mol
        id = names.unique_name()
        # Append mol to named and unnamed stores
        self.mols[kind].append(mol)
        self.mols_named[kind][id] = mol
        return id

    def nearest_binding_site(self, x):
        """The nearest reachable actin binding site"""
        # Get all candidate actins
        actins = [t.mols["actin"] for t in self.reachable]
        actins = list(itertools.chain(*actins))  # flatten
        # Find the g-actin pair nearest our location
        near = [act.nearest(x) for act in actins]
        distances = np.abs(np.subtract([g.x for g in near], x))
        nearest = near[np.argmin(distances)]
        return nearest
