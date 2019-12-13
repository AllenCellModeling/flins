# encoding: utf-8
"""
What is bound to what.

Let different proteins link to and unlink from each other.
"""


class BindingSite:
    """A link between this (thing) and that (thing)"""

    def __init__(self, parent):
        """Create the binding site

        Parameters
        ----------
        parent : `protein`
            Protein on which this binding site is located
        """
        self.parent = parent
        self.link = None

    def __str__(self):
        """String representation of a binding site"""
        if not self.bound:
            return "Unbound site"
        else:
            return "Site bound to " + str(self.linked)

    @property
    def linked(self):
        """Return the parent object of the linked site (None if unbound)"""
        if self.bound:
            return self.link.parent
        else:
            return None

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
