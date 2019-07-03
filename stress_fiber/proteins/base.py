# encoding: utf-8
"""
Base protein class
"""

import numpy as np

from ..base import Base


class Protein(Base):
    def __init__(self, kind, tract=None):
        self.kind = kind
        if tract is not None:
            self._link_tract(tract)
        else:
            self.tract = None
            self.id = None
            self.address = ((self.kind, None),)

    def _link_tract(self, tract):
        """Link a tract, simultaneously creating the local ID and address"""
        self.tract = tract
        self.id = tract.add_mol(self.kind, self)
        self.address = (tract.address[:], (self.kind, self.id))

    @property
    def _space_limits(self):
        """What are the X limits available to this protein?

        Defined as 0 to tract span (if given tract) else 0 to infinity.
        """
        if self.tract is None:
            return (0, np.inf)
        else:
            return (0, self.tract.space.span)
