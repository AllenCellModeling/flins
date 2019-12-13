#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test tracts.
"""

import pytest

import flins.space as space

space_list = [space.Space("hex", r, 100) for r in (0, 1, 2, 3)]
tract_list = [t for s in space_list for t in s.all_tracts]
solo_tract = space.Tract((1, -1, 0), None)


class TestTract:
    @pytest.mark.parametrize("tract", tract_list)
    def test__str__(self, tract):
        """Should not return a <blah>"""
        assert not str(tract).startswith("<")

    @pytest.mark.parametrize("tract", tract_list)
    def test_neighbors(self, tract):
        """Should return a set of neighbors, all of which are distance 1"""
        dist_fn = tract.space.grid.distance
        for neighbor in tract.neighbors:
            dist = dist_fn(tract.loc, neighbor.loc)
            assert dist == 1

    def test_solo_neighbors(self):
        """Should return None"""
        assert solo_tract.neighbors is None

    @pytest.mark.parametrize("tract", tract_list)
    def test_reachable(self, tract):
        """Should return a set of neighbors and self"""
        dist_fn = tract.space.grid.distance
        for reachable in tract.reachable:
            dist = dist_fn(tract.loc, reachable.loc)
            assert dist == 1 or dist == 0

    def test_solo_reachable(self):
        """Should return list of self"""
        assert solo_tract.reachable == [solo_tract]
