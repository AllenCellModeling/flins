#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Test location management
"""

import pytest

import numpy as np
import scipy as sp
import scipy.stats
import stress_fiber as sf
from stress_fiber.construct import locations

np.random.seed(0)
distributions = (
    locations.normal_length_distribution(1000, 10, 10),
    locations.normal_length_distribution(1000, 100, 1),
    locations.normal_length_distribution(1000, 100, 10),
    locations.normal_length_distribution(1000, 200, 50),
    locations.normal_length_distribution(1000, 300, 50),
)


def _variation(locations, lengths, span):
    img = np.zeros((lengths.size, span))
    for i, (length, loc) in enumerate(zip(lengths, locations)):
        img[i, int(loc) : int(loc + length)] = 1
    density = img.sum(0)
    variation = np.std(density) / np.mean(density)
    return variation


@pytest.mark.parametrize("dist", distributions)
def test_normal_length_distribution(dist):
    assert not any(dist < 0), "Locs must be positive"
    with pytest.raises(Exception):
        locations.normal_length_distribution(10, -1, 0)


@pytest.mark.parametrize("dist", distributions)
def test_location_end_pushed(dist):
    span = 3000
    locs, lens = locations.location_end_pushed(dist, span)
    assert all(lens == dist), "Push doesn't mod lengths"


@pytest.mark.parametrize("dist", distributions)
def test_location_end_cut(dist):
    span = 3000
    cut_var = _variation(*locations.location_end_cut(dist, span), span)
    push_var = _variation(*locations.location_end_pushed(dist, span), span)
    assert cut_var < push_var, "Cut has less variation"
    # No good test to prove even density that doesn't fail occasionally...
