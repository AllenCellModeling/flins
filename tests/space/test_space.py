#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test space.
"""

import pytest
import random

import stress_fiber.space as space


space_list = [space.Space("hex", r, 100, True) for r in (0, 1, 2, 3)]
space_list += [space.Space("hex", r, 100, False) for r in (1, 3)]


class TestSpace:
    @pytest.mark.parametrize("space", space_list)
    def test__str__(self, space):
        """Should not return a <blah>"""
        assert not str(space).startswith("<")

    @pytest.mark.parametrize("space", space_list)
    def test_all_tracts(self, space):
        """Tracts equal to the entries in the grid"""
        tracts = space.all_tracts
        grid_entries = [entry["tract"] for entry in space.grid.all_entries]
        assert set(tracts) == set(grid_entries)

    @pytest.mark.parametrize("space", space_list)
    def test_tract(self, space):
        """Return tracts matching our coordinates"""
        assert all([space.tract(tract.loc) == tract for tract in space.all_tracts])

    @pytest.mark.parametrize("space", space_list)
    def test_neighbors(self, space):
        """All neighbors are one away. Only testing one to limit test runtime."""
        tract = random.choice(space.all_tracts)
        loc = tract.loc
        neighbors = space.neighbors(loc)
        assert all([space.grid.distance(n.loc, loc) == 1 for n in neighbors])
