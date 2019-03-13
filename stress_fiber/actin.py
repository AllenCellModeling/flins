# encoding: utf-8
"""
Actin' up
CDW 2019
"""

import numpy as np

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
        self.xlinker = None

    @property
    def x(self):
        """Where are you at? Referenced from parent actin."""
        return self.filament.sites_x[self.index]

    @property
    def bound(self):
        """Are you bound?"""
        return self.xlinker is not None


class Actin:
    """A 1D actin"""

    def __init__(self, x, n, t):
        """An actin at x with n pairs of g-actin on tract t"""
        self.x = x
        self.pairs = n
        # Calculate actin rise and run, store for later use
        # Numbers derived from Howard (2001), Pg 125
        mon_per_poly = 26  # number of g-actin in a thin filament section
        poly_base_length = 72.0  # length of thin filament section in nm
        poly_base_turns = 12.0  # rotations, given two start
        rev = 2 * np.pi  # a single revolution
        self._pitch = poly_base_turns * rev / mon_per_poly  # rad/actin pair
        self._rise = poly_base_length / mon_per_poly  # nm / actin pair
        self._radius = 3  # nm, Howard (2001), Pg 121
        # Create the binding sites
        self.sites_x = [x + i * self._rise for i in range(n)]
        self.sites = [BindingSite(self, index) for index in range(n)]

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

    @property
    def force(self):
        """Let's assume that the filament is really stiff, like really 
        stiff, and then find the total force on it.
        """
        sites = filter(lambda s: s.xlinker is not None, self.sites)
        force = np.sum([site.xlinker.force() for site in sites])
        return force

    def diffuse(self):
        """Move around a bit, see AlphaActinin.diffuse for more explanation"""
        L, r = self.length, self._radius
        f_drag = _diffuse.Drag.Cylinder.long_axis_translation(L, r)
        d_x = _diffuse.Dx(f_drag)
        self.x += d_x
        return d_x

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
