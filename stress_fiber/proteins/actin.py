# encoding: utf-8
"""
Actin' up.

These actin filaments are the backbones of the stress-fiber. They transmit
force, are the common binding target for most other proteins, and move in
response to gradual remodeling. 
"""

import warnings
import numpy as np
import scipy.optimize

from .base import Protein
from ..base import Base
from ..support import units
from ..support import diffuse
from ..support import binding_site


class GActinPair(Base):
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
        self.address = (filament.address[:], ("gactinpair", index))
        self.bs = binding_site.BindingSite(self)

    def __str__(self):
        """String representation of a pair"""
        return "GActin %i with %s" % (self.index, str(self.bs))

    @property
    def x(self):
        """Where are you at? Referenced from parent actin."""
        return self.filament.pairs_x[self.index]

    @property
    def polarity(self):
        """Is the plus end to the right? NB: Myosin is plus-end directed.
        For a review of polarity in the sarcomere, see `this book section`_.

        .. _this book section: https://www.ncbi.nlm.nih.gov/books/NBK9961/
        """
        return self.filament.polarity

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


class Actin(Protein):
    """A 1D actin filament that has binding sites, diffusion behavior, etc."""

    def __init__(self, x, tract=None, n_pair=None, length=None, polarity=None):
        """An actin at x with n pairs of g-actin in a tract.
        
        Parameters
        ----------
        x : `float`
            X location of the actin within the tract. This is where the first
            pair will be located, with all others calculated in reference to it. 
        tract : `stress_fiber.space.Tract`, optional
            1D tract that the actin lives in. Is the parent of the filament in
            organizational hierarchy. 
        n_pair : `int`, optional
            Number of g-actin pairs in the filament. We use this to set filament
            length. Takes precedence over length.
        length : `float`, optional
            Alternate method of setting number of g-actin pairs. Will set
            n_pair to number that results in filament extent closest to length.
        polarity : `boolean`, optional 
            Is the plus end to the right? Governs myosin binding shortening. If
            not passed, we answer None and let the motors figure it out. 
        """
        # Store tract if given, else create placeholders
        super().__init__("actin", tract)
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
            n = round(length / self._rise)
        else:
            raise Exception("must pass n_pair or length to actin on creation")
        # Store locations
        self.n_pairs = n
        self.x = x
        self.pairs_x = self._calc_pairs_x()  # redundant, but here for reminder
        # Create g-actin pairs
        self.pairs = [GActinPair(self, index) for index in range(n)]
        # Store polarity
        assert polarity in (None, True, False), "Polarity is boolean or none"
        self.polarity = polarity

    def __str__(self):
        """String representation of actin"""
        n, l = self.n_pairs, self.length
        xmin, xmax = self.boundaries
        bound = sum([pair.bs.bound for pair in self.pairs])
        unbound = len(self.pairs) - bound
        force, energy = self.force, self.energy
        state = (
            "Actin w/ %i pairs (%.1fnm), between x=%.1f-%.1f, with %i/%i bound/unbound pairs. Force/energy is %.1fpN/%.1fpN*nm"
            % (n, l, xmin, xmax, bound, unbound, force, energy)
        )
        return state

    def _calc_pairs_x(self, start_x=None):
        """Where would our g-actin pairs be for a given starting location?"""
        if start_x is None:
            start_x = self.x
        try:
            return self.__unshifted_locations + start_x
        except AttributeError:
            rise, n = self._rise, self.n_pairs
            self.__unshifted_locations = np.array([rise * i for i in range(n)])
            return self.__unshifted_locations + start_x

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
    def _drag(self):
        """What a. The drag coefficient for this protein with correction factor
        for cytoplasmic sub-diffusion
        """
        try:
            return self.__drag
        except AttributeError:
            L, r = self.length, self._radius
            drag = diffuse.Drag.Cylinder.long_axis_translation(L, r)
            cytoplasmic_subdiffusion = units.world.D_cyto_corr
            self.__drag = drag * cytoplasmic_subdiffusion
        return self.__drag

    @property
    def bound(self):
        """Do you have any attached pairs?"""
        return any([pair.bs.bound for pair in self.pairs])

    def nearest(self, x):
        """What is the nearest pair to the given location"""
        if x < self.x:
            return self.pairs[0]
        elif x > self.x + self.length:
            return self.pairs[-1]
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
        f_drag = diffuse.Drag.Cylinder.long_axis_translation(L, r)
        d_x = diffuse.Dx(f_drag)
        self.x += d_x
        if self.tract is not None:  # then derive diffusion limits from tract
            start, end = self.boundaries
            self.x, _ = diffuse.coerce_to_bounds(start, end, self._space_limits)
        return d_x

    def step(self):
        """Take a timestep: move subject to force and diffusion
        As with α-actinin, we take free diffusion to be subject to an
        Einstein-Smoluchowski diffusive processes (as documented in
        `stress_fiber.support.diffuse`).

        When bound we treat the movement of actin as subject to equipartition of
        energy stored in each of the bound α-actinins. We first find the
        x-location where the force from bound α-actinins is balanced,
        representing relaxation into the local lowest-energy state. We then draw
        a Boltzmann distributed energy, use its sign as a directionality
        indicator, and move in the given direction until the energy stored in
        the system has changed by the drawn amount. 
        """
        """Take a timestep subject to force balance and diffusion"""
        # Find force balanced location
        force_x = self._x_with_balanced_forces
        # Find energy at that location and perturbed energy sampled
        base_energy = self._hypothetical_energy(force_x)
        energy_target = np.random.normal(0, 0.5 * units.constants.kT)
        # Don't move if in energy constrained state already?
        if base_energy >= abs(energy_target):
            self.x = force_x  # update x to energy-locked location
            energy_x = force_x
        else:
            energy_x = self._move_to_dissipate_energy(energy_target, force_x)
            self.x = energy_x
        return (force_x, energy_x)

    @property
    def _x_with_balanced_forces(self):
        """Where would have balanced forces?"""
        # If not bound, here
        if not self.bound:
            return self.x
        # X location limits
        x_limits = (self._space_limits[0], self._space_limits[1] - self.length)
        window_x = lambda x: max(x_limits[0], min(x, x_limits[1]))
        # Balance forces, finding local relaxation point
        optim_out = scipy.optimize.least_squares(self._hypothetical_force, self.x)
        if optim_out.success is not True:
            warnings.warn("Unsuccessful force minimization: " + str(optim_out))
        minimal_force_x = window_x(optim_out.x[0])
        return minimal_force_x

    def _energy_difference(self, x, x_i, energy_budget):
        """Diffusion energy usage calculation

        This is the energy difference between target energy to be consumed and
        energy taken to move this far so far

        Parameters
        ----------
        x: float
            Proposed current location of actin
        x_i: float
            What location we are assuming as the start
        energy_budget: float
            Energy to be consumed during movement
        """
        drag_energy = self._drag * abs(x_i - x)
        spring_energy = self._hypothetical_energy(x)
        mismatch = abs(energy_budget) - (drag_energy + spring_energy)
        return mismatch

    def _move_to_dissipate_energy(self, energy_target, starting_x):
        """Move until energy is dissipated

        Move in the direction specified by the sign of the energy target until
        all the energy used by the drag of the filament and the energy taken up
        by bound springs has risen to the target energy budget. 

        Parameters
        ----------
        energy_target: float
            Move until you hit this energy level
        starting_x: float
            Location to take as start point for drag calculations
        """
        # Find direction of movement and spatial limits
        x_limits = (self._space_limits[0], self._space_limits[1] - self.length)
        if energy_target < 0:  # move left if change is negative
            bounds = (x_limits[0], starting_x)
        else:  # move right if change is positive
            bounds = (starting_x, x_limits[1])
        if bounds[0] == bounds[1]:  # you are pegged at an edge
            return starting_x
        # Record starting energy and solve for new location
        base_energy = self._hypothetical_energy(starting_x)
        energy_least_sq = scipy.optimize.least_squares(
            self._energy_difference,
            starting_x,
            bounds=bounds,
            args=(starting_x, energy_target),
            ftol=np.finfo(float).eps,
        )
        if energy_least_sq.success is not True:
            warnings.warn("Unsuccessful energy minimization: " + str(energy_least_sq))
        minimal_energy_x = energy_least_sq.x[0]
        return minimal_energy_x
