# encoding: utf-8
"""
Cubic coordinate helper functions assuming pointy side up.

This draws heavily from the excellent reference:
    https://www.redblobgames.com/grids/hexagons/
"""

import itertools
import numpy as np

from . import offset
from . import axial


def to_axial(i, j, k):
    """Convert cube coordinates to axial coordinates"""
    q, r = i, k
    return q, r


def to_offset(i, j, k):
    """Convert cube coordinates to offset coordinates, odd rows shifted by +1/2 col"""
    col = i + (j - (j & 1)) // 2
    row = j
    return col, row


def to_cart(i, j, k):
    """Convert cube coordinates to cartesian"""
    # col, row = to_offset(i, j, k)
    # x, y = offset.to_cart(col, row)
    x, y = axial.to_cart(*to_axial(i, j, k))
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
        within_radius = distance(i, j, k, *original) <= n
    return within_radius


def validate(i, j, k):
    """Is this a valid cube coordinate?"""
    valid = (i + j + k) == 0
    return valid


def neighbors(i, j, k):
    """Give me a list of neighboring hexagon locations"""
    neighbors = [
        (i, j + 1, k - 1),
        (i + 1, j, k - 1),
        (i + 1, j - 1, k),
        (i, j - 1, k + 1),
        (i - 1, j, k + 1),
        (i - 1, j + 1, k),
    ]
    return neighbors


def rotate_about_center(i, j, k, n_steps):
    """Rotate the coordinate i,j,k about the 0,0,0 center by n_steps

    Each step is one 60 degree rotation to the right. Negative steps are 
    60 degree rotations to the left.
    """
    coord = np.roll([i, j, k], n_steps)
    if n_steps % 2 == 1:  # odd rotations invert signs
        coord *= -1
    return list(coord)


def closest(i, j, k, points):
    """Which point (in a list) is closest to a single passed point?"""
    dist = lambda pt: distance(i, j, k, *pt)
    dists = [dist(pt) for pt in points]
    return points[np.argmin(dists)]
