# encoding: utf-8
"""
Space is the place we store stuff. 

A space accounts for the interactions of tracts that are arranged according 
to some grid layout. 
"""

from .grids import HexGrid, RectGrid
from .tract import Tract
from ..base import Base


class Space(Base):
    """A spatial region of filaments/crosslinkers

    Each Space is connected to neighbors in a hexagonal grid that gives
    many of the properties of 3D space but with easier mechanics. 
    """

    def __init__(self, kind, size, span=None, mirror=True):
        """ A generic template for 
        Parameters
        ----------
        kind: "hex" or "rect"
            overall shape of the space
        size: int or tuple
            Radius or x/y size of space
        span: None or float
            Length of the tracts in x dimension, optional
        tracts: None or list
            List of tracts
        """
        self.kind = kind
        self.size = size
        self.span = span
        self.mirror = mirror
        if kind == "hex":
            self.grid = HexGrid(size)
        elif kind == "rect":
            self.grid = RectGrid(size)
        else:
            raise Exception("unrecognized tract kind")
        self._tracts = []
        for entry in self.grid.all_entries:
            tract = Tract(entry["cube"], self)
            entry["tract"] = tract
            self._tracts.append(tract)

    def __str__(self):
        """String representation of space"""
        size, n = str(self.size), len(self.all_tracts)
        sp = self.span if self.span is not None else 0
        return "Space with size %s, %i tracts, and span %i" % (size, n, sp)

    @property
    def all_tracts(self):
        """Give all the tracts back"""
        return self._tracts

    def tract(self, loc):
        """Single tract at a given location
        
        Parameters
        ----------
        loc: tuple
            coordinates for the tract we want
        mirror: boolean (True)
            Whether or not to mirror across tractspace edges

        Returns
        -------
        tract or None
        """
        within = self.grid.within(loc)
        valid = self.grid.validate(loc)
        if not valid:
            return None
        if not within:
            if self.mirror:
                loc = self.grid.mirror(loc)
            else:
                return None
        return self.grid.entry(loc)["tract"]

    def neighbors(self, loc):
        """Give me the neighbors to the specified location"""
        ## Validate locations and space size
        within = self.grid.within(loc)
        valid = self.grid.validate(loc)
        if not valid or not within:
            return None
        if self.size == 0:
            return []  # space with only one entry
        ## Find and return neighbors
        n_locs = self.grid.neighbors(loc)
        n_tracts = [self.grid.entry(loc)["tract"] for loc in n_locs]
        n_tracts = list(set(n_tracts))  # dedupe
        return n_tracts
