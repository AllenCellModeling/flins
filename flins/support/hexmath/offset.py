# encoding: utf-8
"""
Grid-scale computing. Offset coordinates.

Hexagonal offset coordinate helper functions, assuming pointy side up.
"""

import numpy as np


def to_cart(col, row):
    """Convert offset coordinates (odd row shifted +1/2 col), into Cartesian
    coordinates.
    NOTE: This produces an output rotated by 60deg from that of `axial.to_cart`
    """
    x = np.sqrt(3) * (col + 0.5 * (row & 1))
    y = 3 / 2 * row
    return x, y


def to_cube(col, row):
    """Convert offset coordinates to cube coordinates"""
    i = col - (row - (row & 1)) // 2
    j = row
    k = -i - j
    return (i, j, k)
