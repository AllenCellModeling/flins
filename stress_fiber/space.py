# encoding: utf-8
"""
Spacing it out
CDW 2019
"""

import itertools
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches


class HexMath:
    """Hexagonal coordinate helper functions, assuming pointy side up.
    This draws heavily from the excellent reference:
        https://www.redblobgames.com/grids/hexagons/
    """

    def cube_to_axial(i, j, k):
        """Convert cube coordinates to axial coordinates"""
        q, r = i, k
        return q, r

    def axial_to_cube(q, r):
        """Convert axial coordinates to cube coordinates"""
        i, k = q, r
        j = -i - k
        return i, j, k

    def cube_to_offset(i, j, k):
        """Convert cube coordinates to offset coordinates, odd rows shifted by +1/2 col"""
        col = i + (k - (k & 1)) // 2
        row = k
        return col, row

    def offset_to_cart(col, row):
        """Convert offset coordinates (odd row shifted +1/2 col), into cartesian coordinates"""
        x = np.sqrt(3) * (col + 0.5 * (row & 1))
        y = 3 / 2 * row
        return x, y

    def axial_to_cart(q, r):
        """Convert axial coordinates into cartesian coordinates"""
        x = np.sqrt(3) * q + np.sqrt(3) / 2 * r
        y = 3 / 2 * r
        return x, y

    def cube_to_cart(i, j, k):
        """Convert cube coordinates to cartesian"""
        col, row = HexMath.cube_to_offset(i, j, k)
        x, y = HexMath.offset_to_cart(col, row)
        return x, y

    def cube_distance(i_1, j_1, k_1, i_2, j_2, k_2):
        """Distance between two hexagons in cube coordinates"""
        a = np.array((i_1, j_1, k_1))
        b = np.array((i_2, j_2, k_2))
        dist = np.sum(np.abs(np.subtract(a, b))) / 2
        return dist

    def cube_within_radius(i, j, k, n, original=(0, 0, 0)):
        """Is this new location within a radius from an original (default origin)"""
        within_radius = HexMath.cube_distance(i, j, k, *original) <= n
        return within_radius

    def cube_validate(i, j, k):
        """Is this a valid cube coordinate?"""
        valid = (i + j + k) == 0
        return valid

    def cube_neighbors(i, j, k):
        """Give me a list of neighboring hexagons"""
        neighbors = [
            (i, j + 1, k - 1),
            (i + 1, j, k - 1),
            (i + 1, j - 1, k),
            (i, j - 1, k + 1),
            (i - 1, j, k + 1),
            (i - 1, j + 1, k),
        ]
        return neighbors

    def cube_to_array_indices(i, j, k, n):
        """Convert cube coordinates to array location for a grid of radius n"""
        x = i + n
        y = k + n
        return x, y

    def create_grid_array(n):
        """Create a list of hex locations for a grid of radius n"""
        # Create grid
        scan = list(range(-n, n + 1))
        grid = [[[] for q in scan] for r in scan]
        # Create tests to ensure location is within roi and valid
        within_radius = lambda i, j, k: HexMath.cube_within_radius(i, j, k, n)
        belongs = lambda ijk: within_radius(*ijk) and HexMath.cube_validate(*ijk)
        # Populate grid with coordinates
        for q, r in itertools.product(scan, scan):
            i, j, k = HexMath.axial_to_cube(q, r)
            if belongs((i, j, k)):
                x, y = HexMath.cube_to_array_indices(i, j, k, n)
                grid[x][y] = {"cube": (i, j, k)}
        grid = np.array(grid)
        return grid


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
