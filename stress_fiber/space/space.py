# encoding: utf-8
"""
Space is the place we store stuff. 

Tracts are sections of space that extend in a long direction and have adjacent
neighboring tracts that also extend in a long direction.
"""

import uuid
import copy
import itertools
import numpy as np

from .hexmath import Cube


class TractSpace:
    """A spatial region of filaments/crosslinkers
    Each TractSpace is connected to neighbors in a hexagonal array that gives
    many of the aspects of 3D space with easier mechanics. 
    """

    def __init__(self, size, span=None):
        """A hexagonal grid of a given symmetric radial size

        Parameters
        ----------
        size: int
            Radius (0 is a single tract, 1 is 9 tracts, ...) of tracts
        span: None or float
            Length of the tracts in x dimension, optional
        """
        self.size = size
        self.span = span
        hex_grid = Cube.create_grid_array(size)
        for n, row in enumerate(hex_grid):
            for m, location in enumerate(row):
                if not location == []:
                    hex_grid[n][m] = Tract(location["cube"], self)
        self._tracts = hex_grid
        self._mirror_centers = Cube.mirrored_centers(size)

    def __str__(self):
        """String representation of tractspace"""
        state = "TractSpace with radius %i, %i tracts" % (
            self.size,
            len(self.all_tracts),
        )
        return state

    def _OOB(self, x, y, z):
        """Is a location out of bounds for list storage?"""
        n = self.size
        return (x + n) > n or (y + n) > n or (z + n) > n

    def tract(self, i, j, k, mirror=True):
        """A single tract at cube coordinates
        
        Parameters
        ----------
        i,j,k: int
            The cube coordinates for the tract we want
        mirror: boolean (True)
            Whether or not to mirror across tractspace edges
        """
        n = self.size
        within = Cube.within_radius(i, j, k, n)
        valid = Cube.validate(i, j, k)
        if not valid:
            return None
        if not within:
            if mirror:
                i, j, k = Cube.mirror(i, j, k, n, self._mirror_centers)
            else:
                return None
        x, y = Cube.to_array_indices(i, j, k, n)
        tract = self._tracts[x, y]
        return tract

    @property
    def all_tracts(self):
        """Give all the tracts back"""
        all_tracts = [tract for tract in self._tracts.flat if not tract == []]
        return all_tracts

    def neighbors(self, i, j, k):
        """Give me the neighbors, ignoring out of bounds"""
        if not Cube.within_radius(i, j, k, self.size):
            return None  # OOB
        if self.size == 0:
            return list()
        neighboring_coordinates = Cube.neighbors(i, j, k)
        tracts = [
            self.tract(i, j, k)
            for i, j, k in neighboring_coordinates
            if self.tract(i, j, k) is not None
        ]
        tracts = list(set(tracts))  # dedupe
        return tracts


class Tract:
    """A single tract in a TractSpace"""

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
                self._neighbors = self.space.neighbors(*self.loc)
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
        if not kind in self.mols:
            self.mols[kind] = []
        if not kind in self.mols_named:
            self.mols_named[kind] = {}
        # Create id for this mol
        id = uuid.uuid4().hex
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
