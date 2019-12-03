# encoding: utf-8
"""
Make a world

Each world tracks the execution of a run. 
"""

import numpy as np


class World:
    """Keep track of a tract space, simulation time, and metadata"""

    def __init__(self, tractspace, random_state=None, **kwargs):
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
        all_mols = [
            m
            for t in self.tractspace.all_tracts
            for v in t.mols.values()
            for m in list(v)
        ]
        for mol in np.random.permutation(all_mols):
            mol.step()
