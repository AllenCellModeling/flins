# encoding: utf-8
"""
_binding_site.py: Track what is bound to what
CDW 2019
"""

class BindingSite:
    """A link between this (thing) and that (thing)"""

    def __init__(self, other=None):
        self.link = None
    
    @property
    def bound(self):
        """Are you currently bound?"""
        return self.link is not None

    def bind(self, other):
        """Bind to other object"""
        assert not self.bound, "Tried to link an already-bound site"
        assert not other.bound, "Tried to link to an already-bound site"
        self.link = other
        self.link.link = self

    def unbind(self):
        """Unbind from other object"""
        assert self.bound, "Tried to unlink an unbound site"
        assert self.link.bound, "Linked site already unbound? Weird."
        self.link.link = None
        self.link = None
