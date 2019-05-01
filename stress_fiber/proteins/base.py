# encoding: utf-8
"""
Base protein class
"""

from ..base import Base


class Protein(Base):
    def __init__(self, kind, tract=None):
        self.kind = kind
        if tract is not None:
            self._link_tract(tract)
        else:
            self.tract = None
            self.id = None
            self.address = None

    def _link_tract(self, tract):
        """Link a tract, simultaneously creating the local ID and address"""
        self.tract = tract
        self.id = tract.add_mol(self.kind, self)
        self.address = (tract.address, (self.kind, self.id))
