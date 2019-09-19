# encoding: utf-8
"""
Grid-scale computing. 

Hexagonal offset coordinate helper functions, assuming pointy side up.
"""

import itertools
import numpy as np


def to_cart(col, row):
    """Convert offset coordinates (odd row shifted +1/2 col), into Cartesian coordinates"""
    x = np.sqrt(3) * (col + 0.5 * (row & 1))
    y = 3 / 2 * row
    return x, y
