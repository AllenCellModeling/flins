# encoding: utf-8

"""
Functions that support the creation of new worlds
"""

import numpy as np


def normal_length_distribution(n, peak, scale):
    """Normally distributed set of lengths, negative values thrown out

    Parameters
    ----------
    n: int
        Number of lengths to return
    peak: float
        Mean length value (assuming no rejected negative values)
    scale: float
        Spread of length value
    """
    if peak < 0 or scale <= 0:
        raise Exception("Length peak and scale must be positive")
    lengths = np.random.normal(loc=peak, scale=scale, size=n)
    if any(lengths < 0):
        replace = np.count_nonzero(lengths < 0)
        lengths[lengths < 0] = normal_length_distribution(replace, peak, scale)
    return lengths


def location_end_pushed(lengths, span):
    """Randomly distribute along span, forcing end overlaps inwards
    If a protein is poking out of the span (negative or past the span value),
    push it back into the valid range. This maintains the length distribution
    but gives an uneven protein density across space.

    Parameters
    ----------
    lengths: array of floats
        Length of each protein
    span: float
        Span of the space
    """
    locations = (lengths + span) * np.random.random(lengths.size) - lengths
    locations = np.clip(locations, 0, span - lengths)
    return locations, lengths


def location_end_cut(lengths, span):
    """Randomly distribute along span, shorten end proteins if they protrude
    If a protein is poking out of the span, adjust its length so that it no
    longer does. This maintains a constant protein density across space at
    the cost of altering the distribution of lengths. 

    Parameters
    ----------
    lengths: array of floats
        Length of each protein
    span: float
        Span of the space
    """
    locations = (lengths + span) * np.random.random(lengths.size) - lengths
    lengths[locations < 0] += locations[locations < 0]
    ends = locations + lengths
    lengths[ends > span] -= locations[ends > span] + lengths[ends > span] - span
    locations = np.clip(locations, 0, span - lengths)
    return locations, lengths
