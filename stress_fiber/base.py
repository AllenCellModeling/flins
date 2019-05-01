# encoding: utf-8
"""
Base classes we use to create others
"""


class Base:
    """Global base class inherited by most classes above"""

    def __init__(self):
        return

    def _repr_pretty_(self, p, cycle):
        """Give pretty iPython printing"""
        p.text(str(self) if not cycle else "...")
