# encoding: utf-8
"""
Grid-scale computing. 

Support the calculation of adjacency and distance on a hexagonal grid given
various coordinate systems.
"""

import itertools
import numpy as np


class Cube:
    """Hexagonal coordinate helper functions, assuming pointy side up.
    This draws heavily from the excellent reference:
        https://www.redblobgames.com/grids/hexagons/
    """

    def to_axial(i, j, k):
        """Convert cube coordinates to axial coordinates"""
        q, r = i, k
        return q, r

    def to_offset(i, j, k):
        """Convert cube coordinates to offset coordinates, odd rows shifted by +1/2 col"""
        col = i + (k - (k & 1)) // 2
        row = k
        return col, row

    def to_cart(i, j, k):
        """Convert cube coordinates to cartesian"""
        col, row = Cube.to_offset(i, j, k)
        x, y = Offset.to_cart(col, row)
        return x, y

    def to_array_indices(i, j, k, n):
        """Convert cube coordinates to array location for a grid of radius n"""
        x = i + n
        y = k + n
        return x, y

    def distance(i_1, j_1, k_1, i_2, j_2, k_2):
        """Distance between two hexagons in cube coordinates"""
        a = np.array((i_1, j_1, k_1))
        b = np.array((i_2, j_2, k_2))
        dist = np.sum(np.abs(np.subtract(a, b))) / 2
        return dist

    def within_radius(i, j, k, n, original=(0, 0, 0)):
        """Is this new location within a radius from an original (default origin)"""
        if original == (0, 0, 0):
            within_radius = max((abs(i), abs(j), abs(k))) <= n
        else:
            within_radius = Cube.distance(i, j, k, *original) <= n
        return within_radius

    def validate(i, j, k):
        """Is this a valid cube coordinate?"""
        valid = (i + j + k) == 0
        return valid

    def neighbors(i, j, k):
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

    def rotate_about_center(i, j, k, n):
        """Rotate the coordinate i,j,k about the 0,0,0 center by n steps

        Each step is one 60 degree rotation to the right. Negative steps are 
        60 degree rotations to the left.
        """
        coord = np.roll([i, j, k], n)
        if n % 2 == 1:  # odd rotations invert signs
            coord *= -1
        return list(coord)

    def mirrored_centers(r):
        """Return the mirrored centers of a world with radius r

        These are the locations that the centers of worlds of equal radius would
        occupy in a greater-world coordinate system. These are used to calculate
        wrapping when walking off the edge of the world. 
        """
        i, j, k = 2 * r + 1, -r, -r - 1
        mirrored = [Cube.rotate_about_center(i, j, k, n) for n in range(6)]
        return mirrored

    def closest(i, j, k, points):
        """Which point (in a list) is closest to a single passed point?"""
        dist = lambda pt: Cube.distance(i, j, k, *pt)
        dists = [dist(pt) for pt in points]
        return points[np.argmin(dists)]

    def mirror(i, j, k, r, centers=None):
        """Return the location, mirrored across the boundary if needed

        Use pre-calculated mirrored centers if given, else calculate them.
        Reference here: https://gamedev.stackexchange.com/questions/137603
        """
        ## Before all else, if within radius of world no mirroring is needed
        if Cube.within_radius(i, j, k, r):
            return i, j, k
        ## Calculate mirrored centers if not given
        if centers is None:
            centers = Cube.mirrored_centers(r)
        ## Find closest mirrored center
        closest = Cube.closest(i, j, k, centers)
        ## Subtract that center to shift back into world
        mirrored = np.subtract((i, j, k), closest)
        return list(mirrored)

    def mirrored_distance(i_1, j_1, k_1, i_2, j_2, k_2, r, centers=None):
        """Return the shortest distance between points, allowing mirroring
        """
        # Calculate mirrored centers if not given
        if centers is None:
            centers = Cube.mirrored_centers(r)
        # Calculate mirrored versions of pt 2
        versions = [(i_2, j_2, k_2)]
        for i_h, j_h, k_h in centers:
            versions.append((i_h + i_2, j_h + j_2, k_h + k_2))
        # Find closest and distance
        closest = Cube.closest(i_1, j_1, k_1, versions)
        return Cube.distance(i_1, j_1, k_1, *closest)

    def create_grid_array(n):
        """Create a list of hex locations for a grid of radius n"""
        # Create grid
        scan = list(range(-n, n + 1))
        grid = [[[] for q in scan] for r in scan]
        # Create tests to ensure location is within roi and valid
        within_radius = lambda i, j, k: Cube.within_radius(i, j, k, n)
        belongs = lambda ijk: within_radius(*ijk) and Cube.validate(*ijk)
        # Populate grid with coordinates
        for q, r in itertools.product(scan, scan):
            i, j, k = Axial.to_cube(q, r)
            if belongs((i, j, k)):
                x, y = Cube.to_array_indices(i, j, k, n)
                grid[x][y] = {"cube": (i, j, k)}
        grid = np.array(grid)
        return grid


class Axial:
    """Hexagonal axial coordinate helper functions, assuming pointy side up."""

    def to_cube(q, r):
        """Convert axial coordinates to cube coordinates"""
        i, k = q, r
        j = -i - k
        return i, j, k

    def to_cart(q, r):
        """Convert axial coordinates into cartesian coordinates"""
        x = np.sqrt(3) * q + np.sqrt(3) / 2 * r
        y = 3 / 2 * r
        return x, y


class Offset:
    """Hexagonal offset coordinate helper functions, assuming pointy side up."""

    def to_cart(col, row):
        """Convert offset coordinates (odd row shifted +1/2 col), into Cartesian coordinates"""
        x = np.sqrt(3) * (col + 0.5 * (row & 1))
        y = 3 / 2 * row
        return x, y
