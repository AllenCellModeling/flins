# encoding: utf-8
"""
Spring forward, an elastic helper

Provide a linear spring that can be subject to thermal forcing.
"""

from . import units

from numpy import pi, sqrt
import math as m
import numpy.random as random


class Spring:
    """A generic one-state spring"""

    def __init__(self, k, rest):
        """Give me a spring

        The spring can be extensional or torsional, the module works in both
        cases. 

        Parameters
        ----------
        k : float
            Spring constant in pN per nm
        rest : float
            Rest angle or length of spring in radians or nm
        """
        ## Passed variables
        self.k = k
        self.rest = rest
        ## Diffusion governors
        kT = units.constants.boltzmann * units.constants.temperature
        # Normalize: a factor used to normalize the PDF of the segment values
        self._normalize = sqrt(2 * pi * kT / self.k)
        self._stand_dev = sqrt(kT / self.k)  # of segment values

    def energy(self, length):
        """Given a current length/angle, return stored energy

        Parameters
        ----------
        length : float
            a spring length or angle in radians or nm

        Returns
        -------
        energy
            the energy required to achieve the given value
        """
        return 0.5 * self.k * m.pow((length - self.rest), 2)

    def force(self, length):
        """Given a current length/angle, return the force/torque.

        This is given signs as though the spring is attached at origin 
        and we are looking at the force needed to hold the spring in place.

        Parameters
        ----------
        length : float
            a spring length or angle in radians or nm

        Returns
        -------
        force
            force exerted by the spring at current length/angle
        """
        return self.k * (length - self.rest)

    def bop_dx(self):
        """Bop for a displacement from rest length
        Assume an exponential energy distribution.
        """
        return random.normal(0, self._stand_dev)

    def bop_length(self):
        """Bop for a new spring length that differs from rest value by dx
        Assume an exponential energy distribution.
        """
        return self.rest + self.bop_dx()
