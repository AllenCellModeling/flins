# encoding: utf-8
"""
Make a world

A whole new world in which we track the execution of a model run.
"""

import itertools
import numpy as np

from . import space
from . import proteins


def create_test_world(radius, span, n_actin, n_actinin):
    """Create a world of given radius with n_actin and n_actinin per tract"""
    tractspace = space.TractSpace(radius, span)
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
                proteins.Anchor(x_end, actin.pairs[-1])
        for _ in range(n_actinin):
            x = np.random.rand() * (span - 35)
            proteins.AlphaActinin(x, tract)
    world = World(tractspace)
    return world


class World:
    """Keep track of a tract space, simulation time, and metadata"""

    def __init__(self, tractspace, random_state=None):
        """Save the tractspace, initiate record keeping 
        
        Parameters
        ----------
        tractspace : `stress_fiber.space.space.TrackSpace`
            Spatial component of the world
        random_state : 
            The internal state of numpy's Mersenne Twister implementation as
            given by `numpy.random.get_state()`. This allows us to recreate run
            trajectories. 
        """
        if random_state is None:
            np.random.seed()  # Ensure proper seeding, making each world unique
            self._starting_random_state = np.random.get_state()
        else:
            np.random.set_state(random_state)
            self._starting_random_state = random_state
        self.tractspace = tractspace
        self.time = 0

    def step(self):
        """Step forward one tick"""
        self.time += 1
        for tract in np.random.permutation(self.tractspace.all_tracts):
            all_mols = list(itertools.chain(*tract.mols.values()))
            for mol in np.random.permutation(all_mols):
                mol.step()
