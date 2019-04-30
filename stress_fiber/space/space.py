# encoding: utf-8
"""
Space is the place we store stuff. 

Tracts are sections of space that extend in a long direction and have adjacent
neighboring tracts that also extend in a long direction.
"""

from .hexmath import HexMath

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches


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
        hex_grid = HexMath.create_grid_array(size)
        for n, row in enumerate(hex_grid):
            for m, location in enumerate(row):
                if not location == []:
                    hex_grid[n][m] = Tract(location["cube"], self)
        self._tracts = hex_grid

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

    def tract(self, i, j, k):
        """A single tract at cube coordinates"""
        n = self.size
        within = HexMath.cube_within_radius(i, j, k, n)
        valid = HexMath.cube_validate(i, j, k)
        if not within or not valid:
            return None  # OOB
        x, y = HexMath.cube_to_array_indices(i, j, k, n)
        tract = self._tracts[x, y]
        return tract

    @property
    def all_tracts(self):
        """Give all the tracts back"""
        all_tracts = [tract for tract in self._tracts.flat if not tract == []]
        return all_tracts

    def neighbors(self, i, j, k):
        """Give me the neighbors, ignoring out of bounds"""
        if not HexMath.cube_within_radius(i, j, k, self.size):
            return None  # OOB
        neighboring_coordinates = HexMath.cube_neighbors(i, j, k)
        tracts = [
            self.tract(i, j, k)
            for i, j, k in neighboring_coordinates
            if self.tract(i, j, k) is not None
        ]
        return tracts

    def plot(self, callback=None):
        """Plot the tracts and what's in them, optionally using a callback"""
        # Set up callback
        if callback is None:
            callback = lambda t: t["cube"]  # default to cube coordinates
        # Set up figure
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.axis("off")
        limit = 2 * self.size + 1
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect=1)
        # Work through each tract
        for tract in self.all_tracts:
            x, y = HexMath.cube_to_cart(*tract["cube"])
            hex = matplotlib.patches.RegularPolygon(
                (x, y), numVertices=6, radius=1, facecolor="White", edgecolor="k"
            )
            ax.add_patch(hex)
            ax.text(x, y, str(callback(tract)), ha="center", va="center")
        plt.show()


class Tract:
    """A single tract in a TractSpace"""

    def __init__(self, loc, space):
        """A single tract in a space"""
        self.loc = loc
        self.space = space
        self._neighbors = None
        self.mols = {}

    def __str__(self):
        """String representation of tract"""
        loc_str = str(self.loc)
        mols_counts = ["%i %ss" % (len(v), k) for k, v in self.mols.items()]
        mols_str = ", ".join(mols_counts)
        return "Tract at %s with %s" % (loc_str, mols_str)

    @property
    def neighbors(self):
        """Who are your neighbors?
        We don't want to calculate this each time, but we can't populate it on
        creation as there are neighbors that haven't yet been created. So we
        populate this list on the first call and then reference the stored
        version thereafter."""
        if self._neighbors is None:
            self._neighbors = self.space.neighbors(*self.loc)
        return self._neighbors
