#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Test α-actinin
"""

import pytest

import numpy as np
import stress_fiber as sf
from stress_fiber.proteins.alpha_actinin import AlphaActinin


def actinin_tract():
    span = 10000
    tractspace = sf.space.Space("hex", 0, span)
    t = tractspace.all_tracts[0]
    t.rand = lambda: span * (np.random.rand() * 0.5 + 0.25)
    return t


def unbound_actinin():
    t = actinin_tract()
    _ = sf.proteins.Actin(t.rand(), t, length=100)
    return AlphaActinin(t.rand(), t)


def one_side_bound_actinin(head_i=0):
    t = actinin_tract()
    alph = AlphaActinin(t.rand(), t)
    act = sf.proteins.Actin(alph.heads[head_i].x, t, length=100)
    gact = act.pairs[0]
    alph.heads[head_i].bs.bind(gact.bs)
    return alph


def both_sides_bound_actinin():
    t = actinin_tract()
    alph = AlphaActinin(t.rand(), t)
    acts = [sf.proteins.Actin(head.x, t, length=100) for head in alph.heads]
    gacts = acts[0].pairs[0], acts[1].pairs[0]
    alph.heads[0].bs.bind(gacts[0].bs)
    alph.heads[1].bs.bind(gacts[1].bs)
    return alph


actinin_list = [
    unbound_actinin(),
    both_sides_bound_actinin(),
    one_side_bound_actinin(0),
    one_side_bound_actinin(1),
]

head_list = [head for actinin in actinin_list for head in actinin.heads]


# Testing begins


def test_actinin_world():
    """Create a world with some actinins and run it for a bit"""
    w = sf.construct.create_test_world(1, 1000, 2, 10, 0)
    _ = [w.step() for i in range(10)]
    tracts = w.tractspace.all_tracts
    actinins = [a for t in tracts for a in t.mols["actinin"]]
    states = [a.bound for a in actinins]
    assert not all(states) and any(states)  # mix of bound/unbound


class TestActinin:
    @pytest.mark.parametrize("actinin", actinin_list)
    def test__str__(self, actinin):
        """Are our nice string representations getting called?"""
        assert not str(actinin).startswith("<")

    @pytest.mark.parametrize("actinin", actinin_list)
    def test_bound(self, actinin):
        """Are our actinins bound?"""
        assert actinin.bound == any([h.bs.bound for h in actinin.heads])

    @pytest.mark.parametrize("actinin", actinin_list)
    def test_fully_bound(self, actinin):
        """Are our actinins fully bound?"""
        assert actinin.fully_bound == all([h.bs.bound for h in actinin.heads])

    @pytest.mark.parametrize("actinin", actinin_list)
    def test_energy(self, actinin):
        """Are our energies only present when fully bound?"""
        if not actinin.fully_bound:
            assert actinin.energy == 0
        else:
            assert actinin.energy != 0


class TestActininHead:
    @pytest.mark.parametrize("head", head_list)
    def test__str__(self, head):
        """Are our nice string representations getting called?"""
        assert not str(head).startswith("<")

    @pytest.mark.parametrize("head", head_list)
    def test_other_head(self, head):
        """Other head gives diff head on same protein"""
        other_head = head.other_head
        assert head != other_head
        assert head.parent == other_head.parent

    @pytest.mark.parametrize("head", head_list)
    def test__r01(self, head):
        """Ensure that our binding rates are going in the right direction"""
        assert head._r01(-1) < head._r01(0)
        assert head._r01(0) > head._r01(1)
