#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Test α-actinin
"""

import pytest

import numpy as np
import stress_fiber as sf
import stress_fiber.proteins.alpha_actinin as aact


class TestActininHead:
    TRACTSPACE = sf.space.TractSpace(1, 100)
    TRACTS = TRACTSPACE.all_tracts
    ACTININS = [aact.AlphaActinin(np.random.randint(100), t) for t in TRACTS]
    HEADS = [head for actinin in ACTININS for head in actinin.heads]

    @pytest.mark.parametrize('head', HEADS)
    def test__str__(self, head):
        """Are our nice string representations getting called?"""
        assert not str(head).startswith("<")
    
    @pytest.mark.parametrize('head', HEADS)
    def test_other_head(self, head):
        """Other head gives diff head on same protein"""
        other_head = head.other_head
        assert head != other_head
        assert head.actinin == other_head.actinin

    @pytest.mark.parametrize('head', HEADS)
    def test__r01(self, head):
        """Ensure that our binding rates are going in the right direction"""
        assert head._r01(-1) < head._r01(0)
        assert head._r01(0) > head._r01(1)




