# encoding: utf-8
"""
_spring: Spring forward, a helper function
CDW 2019
"""

from . import _units

from numpy import pi, sqrt
import math as m
import numpy.random as random

random.seed()  # Ensure proper seeding


class Spring:
    """A generic one-state spring"""

    def __init__(self, k, rest):
        """Give me a spring"""
        ## Passed variables
        self.k = k
        self.rest = rest
        ## Diffusion governors
        kT = _units.constants.boltzmann * _units.constants.temperature
        # Normalize: a factor used to normalize the PDF of the segment values
        self._normalize = sqrt(2 * pi * kT / self.k)
        self._stand_dev = sqrt(kT / self.k)  # of segment values

    def energy(self, length):
        """Given a current length/angle, return stored energy

        Takes:
            x: a spring length or angle
        Returns:
            energy: the energy required to achieve the given value
        """
        return 0.5 * self.k * m.pow((length - self.rest), 2)

    def force(self, length):
        """Given a current length/angle, return the force/torque.
        This is given signs as though the spring is attached at origin 
        and we are looking at the force needed to hold the spring in place.
        """
        return self.k * (length - self.rest)

    def bop_dx(self):
        """Bop for a displacement from rest length
        Assume an exponential energy dist
        """
        return random.normal(0, self._stand_dev)

    def bop_length(self):
        """Bop for a new spring length that differs from rest length by dx
        Assume an exponential energy dist
        """
        return self.rest + self.bop_dx()
