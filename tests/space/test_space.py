#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test the tract space and tracts. 
"""

import stress_fiber.space as space


class TestTract:
    TRACTSPACES = [space.TractSpace(r, 100) for r in (0, 1, 2, 3)]
    TRACTS = [t for s in TRACTSPACES for t in s.all_tracts]
    SOLO_TRACT = space.Tract((1, -1, 0), None)

    def test_neighbors(self):
        """Should return a set of neighbors, all of which are distance 1"""
        dist_fn = space.hexmath.Cube.mirrored_distance
        for tract in self.TRACTS:
            r = tract.space.size
            for neighbor in tract.neighbors:
                dist = dist_fn(*tract.loc, *neighbor.loc, r)
                assert dist == 1

    def test_solo_neighbors(self):
        """Should return None"""
        assert self.SOLO_TRACT.neighbors is None

    def test_reachable(self):
        """Should return a set of neighbors and self"""
        dist_fn = space.hexmath.Cube.mirrored_distance
        for tract in self.TRACTS:
            r = tract.space.size
            for reachable in tract.reachable:
                dist = dist_fn(*tract.loc, *reachable.loc, r)
                assert dist == 1 or dist == 0

    def test_solo_reachable(self):
        """Should return list of self"""
        assert self.SOLO_TRACT.reachable == [self.SOLO_TRACT]
