# encoding: utf-8
"""
Actin' up
CDW 2019
"""

import warnings
import numpy as np
import scipy.optimize

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches

from . import _units
from . import _diffuse


class BindingSite:
    """An actin binding site"""

    def __init__(self, filament, index):
        self.filament = filament
        self.index = index
        self.link = None

    @property
    def x(self):
        """Where are you at? Referenced from parent actin."""
        return self.filament.sites_x[self.index]

    @property
    def bound(self):
        """Are you bound?"""
        return self.link is not None


class Actin:
    """A 1D actin"""

    def __init__(self, x, n, t):
        """An actin at x with n pairs of g-actin on tract t"""
        # Calculate actin rise and run, store for later use
        # Numbers derived from Howard (2001), Pg 125
        mon_per_poly = 26  # number of g-actin in a thin filament section
        poly_base_length = 72.0  # length of thin filament section in nm
        poly_base_turns = 12.0  # rotations, given two start
        rev = 2 * np.pi  # a single revolution
        self._pitch = poly_base_turns * rev / mon_per_poly  # rad/actin pair
        self._rise = poly_base_length / mon_per_poly  # nm / actin pair
        self._radius = 3  # nm, Howard (2001), Pg 121
        # Store locations 
        self.pairs = n
        self.x = x
        self.sites_x = self._calc_sites_x()  # redundant, but here for reminder
        # Create binding sites
        self.sites = [BindingSite(self, index) for index in range(n)]

    def _calc_sites_x(self, start_x=None):
        """Where would our binding sites be for a given starting location?"""
        if start_x is None:
            start_x = self.x
        rise, n = self._rise, self.pairs
        return np.array([start_x + rise * i for i in range(n)])

    @property
    def x(self):
        """What is the starting location of the actin filament? 
        It is assumed to go in the positive direction afterwards.
        """
        return self._x

    @x.setter
    def x(self, start_x):
        self._x = start_x
        self.sites_x = self._calc_sites_x(start_x)

    @property
    def length(self):
        """How long are you?"""
        return abs(self.pairs * self._rise)

    @property
    def boundaries(self):
        """How far do you extend?"""
        return (self.x, self.x + self.length)

    @property
    def bound(self):
        """Do you have any attached binding sites?"""
        return any([bs.bound for bs in self.sites])

    def nearest(self, x):
        """What is the nearest site to the given location"""
        return self.sites[np.argmin(np.abs(np.subtract(x, self.sites_x)))]

    def _hypothetical_force(self, x):
        """Assume the filament is, like, really stiff and find the force on it.
        
        This takes the sum of forces exerted by bound α-actinins along the
        filament given a hypothetical actin location, x. This is used to do
        force balances without changing the state of the filament.
        """
        bs_x = self._calc_sites_x(x)
        bs_and_x = [(site, x) for site,x in zip(self.sites, bs_x) if site.bound]
        force = np.sum([site.xlinker.force(x) for site,x in bs_and_x])
        return force
        
    @property
    def force(self):
        """How much force does the filament feel at its current location?
        See `_hypothetical_force` for more detail.
        """
        if not self.bound:
            return 0
        force = self._hypothetical_force(self.x)
        #sites = filter(lambda s: s.xlinker is not None, self.sites)
        #force = np.sum([site.xlinker.force() for site in sites])
        return force

    def _hypothetical_energy(self, x):
        """Assume our (axially) stiff actin is storing energy in α-actinin
        What is the current energy in the bound α-actinins?
        """
        bs_x = self._calc_sites_x(x)
        bs_and_x = [(site, x) for site, x in zip(self.sites, bs_x) if site.bound]
        energy = np.sum([site.energy(x) for site, x in bs_and_x])
        return energy

    @property
    def energy(self):
        """How much energy does the filament bear at its current location?"""
        if not self.bound:
            return 0
        energy = self._hypothetical_energy(self.x)
        return energy

    def freely_diffuse(self):
        """Move around a bit, see AlphaActinin.diffuse for more explanation"""
        L, r = self.length, self._radius
        f_drag = _diffuse.Drag.Cylinder.long_axis_translation(L, r)
        d_x = _diffuse.Dx(f_drag)
        self.x += d_x
        return d_x


    def step(self):
        """Take a timestep: move subject to force and diffusion
        As with α-actinin, we take free diffusion to be subject to an
        Einstein-Smoluchowski diffusive processes (as documented in `_diffuse`).

        When bound we treat the movement of actin as subject to equipartition of
        energy stored in each of the bound α-actinins. We first find the
        x-location where the force from bound α-actinins is balanced,
        representing relaxation into the local lowest-energy state. We then draw
        a Boltzmann distributed energy, use its sign as a directionality
        indicator, and move in the given direction until the energy stored in
        the system has changed by the drawn amount. 
        """
        if not self.bound:
            d_x = self.freely_diffuse()
        else:
            # Balance forces, finding local relaxation point
            starting_x = self.x
            force_least_sq = scipy.optimize.least_squares(
                self._hypothetical_force, 
                starting_x)
            minimal_force_x = force_least_sq.x[0]
            if force_least_sq.success is not True:
                warnings.warn("Unsuccessful force minimization: " +
                              str(force_least_sq))
            # Find base and perturbed energies
            base_energy = self._hypothetical_energy(minimal_force_x)
            energy_bump = np.random.normal(0, 0.5 * _units.constants.kT)
            if energy_bump < 0:  # move left if bump is negative
                bounds = (-np.inf, minimal_force_x)
            else:  # move right if bump is positive
                bounds = (minimal_force_x, np.inf)
            energy_bump = abs(energy_bump)
            # Find energy difference and move
            def energy_delta(x):
                energy_from_movement = abs(self._hypothetical_energy(x) - base_energy) 
                energy_mismatch = energy_from_movement - energy_bump
                return energy_mismatch
            energy_least_sq = scipy.optimize.least_squares(
                energy_delta,
                minimal_force_x,
                bounds=bounds)
            if energy_least_sq.success is not True:
                warnings.warn("Unsuccessful energy minimization: " +
                              str(energy_least_sq))
            minimal_energy_x = energy_least_sq.x[0]
            d_x = minimal_energy_x - starting_x
        self.x += d_x

    def plot(self, ax=None, show=False, y=0):
        """Plot this fil"""
        if ax is None:
            fig, ax = plt.subplots(1, 1, figsize=(8, 8))
            ax.axis("off")
            ax.set(
                xlim=self.boundaries,
                ylim=(y - 2 * self._rise, y + 2 * self._rise),
                aspect=1,
            )
        circ = lambda x, y: matplotlib.patches.CirclePolygon(
            (x, y),
            radius=0.66 * self._rise,
            resolution=20,
            facecolor="skyblue",
            edgecolor="royalblue",
        )
        # Work through each site pair
        for x in [s.x for s in self.sites]:
            ax.add_patch(circ(x + self._rise * 0.1, y + 0.5 * self._rise))
            ax.add_patch(circ(x - self._rise * 0.1, y - 0.5 * self._rise))
        if show:
            plt.show()
        return ax
