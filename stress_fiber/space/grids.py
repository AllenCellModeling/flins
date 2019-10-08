# encoding: utf-8
"""
Be a spatial divider, not a uniter

Support the calculation of adjacency, distance, and mirroring on hexagonal 
grids with various shapes.
"""

import abc
from abc import abstractmethod
import itertools
import numpy as np

from ..support.hexmath import cube, axial, offset
from ..base import Base


class Grid(Base, metaclass=abc.ABCMeta):
    @abstractmethod
    def within(self, loc):
        """Is loc within the grid?"""
        pass

    @abstractmethod
    def validate(self, loc):
        """Is loc valid within the grid's coordinate system?"""
        pass

    @abstractmethod
    def mirror(self, loc):
        """Return loc, mirrored across boundaries if outside the grid"""
        pass

    @abstractmethod
    def neighbors(self, loc):
        """Return the coordinates of neighbors"""
        pass

    @abstractmethod
    def distance(self, loc1, loc2):
        """The distance between the locs on this grid"""
        pass

    @abstractmethod
    def entry(self, loc):
        """Give the grid entry at the passed loc"""
        pass

    @property
    @abstractmethod
    def all_entries(self):
        """Return all entries in the grid"""
        pass


class HexGrid(Grid):
    """A world that is in a hexagonal shape with a radius of n hexes"""

    def __init__(self, n, mirror=True):
        """Create a hexagonal grid that we'll use as a model of tracts

        Parameters
        ----------
        n: int
            radius of the world. 0 is a single tract, 1 is seven and so on
        mirror: bool (True)
            whether to mirror when calculating neighbors and distances
        """
        self.size, self._mirroring = n, mirror
        self.array = self._create_grid_array(n)
        if mirror:
            self._mirrored_centers = self._find_mirrored_centers()

    def __str__(self):
        typesize = "HexGrid of radius %i. " % self.size
        mirrored = "Is mirrored." if self._mirroring else "Isn't mirrored."
        return typesize + mirrored

    def within(self, loc):
        """Is the passed loc within the non-mirrored grid?"""
        return cube.within_radius(*loc, self.size)

    def validate(self, loc):
        """Is the passed loc a valid one in cube coords?"""
        return cube.validate(*loc)

    def _find_mirrored_centers(self):
        """Find the mirrored centers of the hexagonal grid

        These are the locations that the centers of worlds of equal radius would
        occupy in a greater-world coordinate system. These are used to calculate
        wrapping when walking off the edge of the world. 
        """
        loc = (2 * self.size + 1, -self.size, -self.size - 1)
        mirrored = [cube.rotate_about_center(*loc, n) for n in range(6)]
        return mirrored

    def mirror(self, loc):
        """Return the location, mirrored across the boundary if needed

        Use pre-calculated mirrored centers if given, else calculate them.
        Reference here: https://gamedev.stackexchange.com/questions/137603
        """
        ## Before all, if within radius of world no mirroring is needed
        if self.within(loc):
            return loc
        ## Find closest mirrored center
        closest = cube.closest(*loc, self._mirrored_centers)
        ## Subtract that center to shift back into world
        mirrored = np.subtract(loc, closest)
        return list(mirrored)

    def neighbors(self, loc: tuple):
        """Return neighboring coordinates, mirroring if desired"""
        if self.size == 0:
            return []
        neighbors = cube.neighbors(*loc)
        assert all([self.validate(n) for n in neighbors])
        if self._mirroring:
            neighbors = [self.mirror(loc) for loc in neighbors]
        else:
            neighbors = [loc for loc in neighbors if self.within(loc)]
        return neighbors

    def distance(self, loc1, loc2):
        """Return the distance between two locations, potentially mirroring"""
        if not self.validate(loc1) or not self.validate(loc2):
            return None
        if not self._mirroring:
            ## Kick out if either are out of bounds
            if not self.within(loc1) or not self.within(loc2):
                return None
        else:
            # Find the mirrored versions of loc2
            centers = self._mirrored_centers
            versions = [loc2]
            for center in centers:
                versions.append(list(np.add(center, loc2)))
            # Find closest and distance
            loc2 = cube.closest(*loc1, versions)
        return cube.distance(*loc1, *loc2)

    def to_array_indices(self, loc):
        """Convert cube coordinates to array location"""
        i, j, k = loc
        x = i + self.size
        y = k + self.size
        return x, y

    def entry(self, loc):
        """Give the grid entry at that location"""
        indices = self.to_array_indices(loc)
        return self.array[indices]

    @property
    def all_entries(self):
        """All entries in the grid"""
        return [entry for row in self.array for entry in row if entry != []]

    def _create_grid_array(self, n):
        """Create a list of hex locations for a grid of radius n"""
        # Create grid
        scan = list(range(-n, n + 1))
        grid = [[[] for q in scan] for r in scan]
        # Create tests to ensure location is within roi and valid
        belongs = lambda loc: self.within(loc) and self.validate(loc)
        # Populate grid with coordinates
        for q, r in itertools.product(scan, scan):
            loc = axial.to_cube(q, r)
            if belongs(loc):
                x, y = self.to_array_indices(loc)
                grid[x][y] = {"cube": loc}
        grid = np.array(grid)
        return grid


class RectGrid(Grid):
    """A rectangular grid of n-by-m hexs"""

    def __init__(self, size, mirror=True):
        """Create a rectangular grid that we'll use as a model of tracts

        Parameters
        ----------
        size: tuple (n, m)
            size of the world, n rows and m columns
        mirror: bool (True)
            whether to mirror when calculating neighbors and distances
        """
        self.size, self._mirroring = size, mirror
        self.array = self._create_grid_array(size)
        if mirror:
            assert size[1] % 2 == 0, "Mirroring only works with even heights"

    def __str__(self):
        typesize = "RectGrid of rows/cols %s. " % str(self.size)
        mirrored = "Is mirrored." if self._mirroring else "Isn't mirrored."
        return typesize + mirrored

    def within(self, loc):
        """Is loc within the non-mirrored grid?"""
        r, c = self.to_array_indices(loc)
        if r >= self.size[0] or c >= self.size[1]:
            return False
        elif r < 0 or c < 0:
            return False
        else:
            return True

    def validate(self, loc):
        """Is loc valid within the grid's coordinate system?"""
        return cube.validate(*loc)

    def mirror(self, loc):
        """Return loc, mirrored across boundaries if outside the grid"""
        rn, cn = self.size
        r, c = self.to_array_indices(loc)
        loc = self._from_array_indices(r % rn, c % cn)
        return loc

    def neighbors(self, loc):
        """Return the coordinates of neighbors"""
        if self.size == (0, 0):
            return []
        neighbors = cube.neighbors(*loc)
        if self._mirroring:
            neighbors = [self.mirror(loc) for loc in neighbors]
        else:
            neighbors = [loc for loc in neighbors if self.within(loc)]
        return neighbors

    def _mirrored_locs(self, loc):
        """Find the real and ghost versions of loc"""
        # Set it up
        rn, cn = self.size
        loc = self.mirror(loc)
        # Doing all this in offset/xy coordinates
        rc = np.array(self.to_array_indices(loc))
        # Wrap around each dimension
        wraps = ((0, 0), (-rn, 0), (rn, 0), (0, -cn), (0, cn), (-rn, -cn), (rn, cn))
        rc_locs = [rc + w for w in wraps]
        locs = [self._from_array_indices(*rc) for rc in rc_locs]
        return locs

    def distance(self, loc1, loc2):
        """The distance between the locs on this grid"""
        if not self.validate(loc1) or not self.validate(loc2):
            return None
        if not self._mirroring:
            ## Kick out if either are out of bounds
            if not self.within(loc1) or not self.within(loc2):
                return None
        else:
            ## Find the wrapped versions of loc2
            versions = self._mirrored_locs(loc2)
            loc2 = cube.closest(*loc1, versions)
        return cube.distance(*loc1, *loc2)

    def entry(self, loc):
        """Give the grid entry at the passed loc"""
        indices = self.to_array_indices(loc)
        return self.array[indices]

    @property
    def all_entries(self):
        """Return all entries in the grid"""
        return [entry for row in self.array for entry in row]

    @staticmethod
    def to_array_indices(loc):
        """Convert to array indices from cube location"""
        i, j, k = loc
        r = k
        # c = -k - j - k // 2
        c = i + k // 2
        return c, r

    @staticmethod
    def _from_array_indices(r, c):
        """Convert array location to cube coordinates"""
        i = r - c // 2
        k = c
        j = -i - k
        return i, j, k

    def _create_grid_array(self, size):
        """Create a grid of hex locations with width, height of size"""
        row_n, col_n = size
        grid = []
        for r in range(row_n):
            cc = range(col_n)
            row = [{"cube": self._from_array_indices(r, c)} for c in cc]
            grid.append(np.array(row))
        return np.array(grid)
