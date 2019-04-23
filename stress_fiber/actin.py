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
from . import _binding_site


class GActinPair:
    """Two g-actin with a site that can bind and unbind."""

    def __init__(self, filament, index):
        """A binding site on actin filament at location index. 

        Parameters
        ----------
        filament : `stress_fiber.Actin`
            Parent actin filament to this pair.
        index : `int`
            Pair location in the list of pairs on the parent filament. 
        """
        self.filament = filament
        self.index = index
        self.bs = _binding_site.BindingSite(self)

    def __str__(self):
        """String representation of a pair"""
        return "GActin %i with %s" % (self.index, str(self.bs))

    @property
    def x(self):
        """Where are you at? Referenced from parent actin."""
        return self.filament.pairs_x[self.index]

    def force(self, x=None):
        """What force does the g-actin experience?"""
        if not self.bs.bound:
            return 0
        if x is None:
            x = self.x
        return self.bs.linked.force(x)

    def energy(self, x=None):
        """What energy does this g-actin store?

        Energy is only stored in the case that the binding site is attached,
        since we are treating actin as inflexible.
        """
        if not self.bs.bound:
            return 0
        if x is None:
            x = self.x
        return self.bs.linked.energy(x)


class Actin:
    """A 1D actin filament that has binding sites, diffusion behavior, etc."""

    def __init__(self, x, tract, n_pair=None, length=None):
        """An actin at x with n pairs of g-actin in a tract.
        
        Parameters
        ----------
        x : `float`
            X location of the actin within the tract. This is where the first
            pair will be located, with all others calculated in reference to it. 
        tract : `stress_fiber.space.tract`
            1D tract that the actin lives in. Is the parent of the filament in
            organizational hierarchy. 
        n_pair : `int`, optional
            Number of g-actin pairs in the filament. We use this to set filament
            length. Takes precedence over length.
        length : `float`, optional
            Alternate method of setting number of g-actin pairs. Will set
            n_pair to number that results in filament extent closest to length.
        """
        # Calculate actin rise and run, store for later use
        # Numbers derived from Howard (2001), Pg 125
        mon_per_poly = 26  # number of g-actin in a thin filament section
        poly_base_length = 72.0  # length of thin filament section in nm
        poly_base_turns = 12.0  # rotations, given two start
        rev = 2 * np.pi  # a single revolution
        self._pitch = poly_base_turns * rev / mon_per_poly  # rad/actin pair
        self._rise = poly_base_length / mon_per_poly  # nm / actin pair
        self._radius = 3  # nm, Howard (2001), Pg 121
        # Calculate length
        if n_pair is not None:
            n = n_pair
        elif length is not None:
            n = round(length/self._rise)
        else:
            raise Exception("must pass n_pair or length to actin on creation")
        # Store locations
        self.n_pairs = n
        self.x = x
        self.pairs_x = self._calc_pairs_x()  # redundant, but here for reminder
        self.tract = tract
        # Create g-actin pairs
        self.pairs = [GActinPair(self, index) for index in range(n)]

    def __str__(self):
        """String representation of actin"""
        n, l = self.n_pairs, self.length
        xmin, xmax = self.boundaries
        bound = sum([pair.bs.bound for pair in self.pairs])
        unbound = len(self.pairs) - bound
        force, energy = self.force, self.energy
        state = (
            "Actin w/ %i pairs (%.1fnm), between x=%.1f-%.1f, \n \
            with %i/%i bound/unbound pairs. Force/energy is %.1fpN/%.1fpN*nm"
            % (n, l, xmin, xmax, bound, unbound, force, energy)
        )
        return state

    def _calc_pairs_x(self, start_x=None):
        """Where would our g-actin pairs be for a given starting location?"""
        if start_x is None:
            start_x = self.x
        rise, n = self._rise, self.n_pairs
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
        self.pairs_x = self._calc_pairs_x(start_x)

    @property
    def length(self):
        """How long are you?"""
        return abs(self.n_pairs * self._rise)

    @property
    def boundaries(self):
        """How far do you extend?"""
        return (self.x, self.x + self.length)

    @property
    def bound(self):
        """Do you have any attached pairs?"""
        return any([pair.bs.bound for pair in self.pairs])

    def nearest(self, x):
        """What is the nearest pair to the given location"""
        return self.pairs[np.argmin(np.abs(np.subtract(x, self.pairs_x)))]

    def nearest_unbound(self, x):
        """What is the nearest unbound pair to the given location"""
        # TODO Implement
        raise NotImplementedError

    def _hypothetical_force(self, x):
        """Assume the filament is, like, really stiff and find the force on it.
        
        This takes the sum of forces exerted by bound α-actinins along the
        filament given a hypothetical actin location, x. This is used to do
        force balances without changing the state of the filament.
        """
        px = self._calc_pairs_x(x)
        pairs_and_x = [(pair, x) for pair, x in zip(self.pairs, px) if pair.bs.bound]
        force = np.sum([pair.force(x) for pair, x in pairs_and_x])
        return force

    @property
    def force(self):
        """How much force does the filament feel at its current location?
        See `_hypothetical_force` for more detail.
        """
        if not self.bound:
            return 0
        force = self._hypothetical_force(self.x)
        return force

    def _hypothetical_energy(self, x):
        """Assume our (axially) stiff actin is storing energy in α-actinin
        What is the current energy in the bound α-actinins?
        """
        px = self._calc_pairs_x(x)
        pairs_and_x = [(pair, x) for pair, x in zip(self.pairs, px) if pair.bs.bound]
        energy = np.sum([pair.energy(x) for pair, x in pairs_and_x])
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
            self.x += d_x
            return
        # Balance forces, finding local relaxation point
        starting_x = self.x
        force_least_sq = scipy.optimize.least_squares(
            self._hypothetical_force, starting_x
        )
        if force_least_sq.success is not True:
            warnings.warn("Unsuccessful force minimization: " + str(force_least_sq))
        minimal_force_x = force_least_sq.x[0]
        # Find base and perturbed energies
        base_energy = self._hypothetical_energy(minimal_force_x)
        energy_bump = np.random.normal(0, 0.5 * _units.constants.kT)
        # Don't move if in energy constrained state already?
        if base_energy >= abs(energy_bump):
            self.x = minimal_force_x
            return (minimal_force_x, minimal_force_x)
        else:
            energy_bump = np.sign(energy_bump) * (base_energy - abs(energy_bump))
        if energy_bump < 0:  # move left if bump is negative
            bounds = (-np.inf, minimal_force_x)
        else:  # move right if bump is positive
            bounds = (minimal_force_x, np.inf)
        energy_bump = abs(energy_bump)
        # Find energy difference and move
        # NOTE: most time intensive part; would benefit from optimization
        def energy_delta(x):
            energy_from_movement = abs(self._hypothetical_energy(x) - base_energy)
            energy_mismatch = energy_from_movement - energy_bump
            return energy_mismatch

        energy_least_sq = scipy.optimize.least_squares(
            energy_delta, minimal_force_x, bounds=bounds
        )
        if energy_least_sq.success is not True:
            warnings.warn("Unsuccessful energy minimization: " + str(energy_least_sq))
        minimal_energy_x = energy_least_sq.x[0]
        self.x = minimal_energy_x
        return (minimal_force_x, minimal_energy_x)

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
        # Work through each pair
        for x in [p.x for p in self.pairs]:
            ax.add_patch(circ(x + self._rise * 0.1, y + 0.5 * self._rise))
            ax.add_patch(circ(x - self._rise * 0.1, y - 0.5 * self._rise))
        if show:
            plt.show()
        return ax
