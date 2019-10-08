# encoding: utf-8
"""
Grid-scale computing. Axial coordinates.

Hexagonal axial coordinate helper functions, assuming pointy side up.
"""

import numpy as np


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
