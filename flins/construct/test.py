# encoding: utf-8
"""
Make a test world

This is the default world used to interactively test model 
components and structure.
"""

import numpy as np

from .world import World
from .. import space
from .. import proteins


def create_test_world(radius, span, n_actin, n_actinin, n_motors):
    """Create a world of given radius with n_actin and n_actinin per tract"""
    tractspace = space.Space("hex", radius, span)
    for tract in tractspace.all_tracts:
        for _ in range(n_actin):
            length = np.random.uniform(0.1 * span, 0.9 * span)
            x = np.random.rand() * (span - length)
            actin = proteins.Actin(x, tract, length=length)
            # Anchor first and last tenth
            if x < span * 0.1:
                proteins.Anchor(x, actin.pairs[0], tract)
            x_end = actin.pairs_x[-1]
            if x_end > span * 0.9:
                proteins.Anchor(x_end, actin.pairs[-1], tract)
        for _ in range(n_actinin):
            x = np.random.rand() * (span - 35)
            proteins.AlphaActinin(x, tract)
        for _ in range(n_motors):
            x = np.random.rand() * (span - 30)
            proteins.Motor(x, tract)
    world = World(tractspace)
    return world
