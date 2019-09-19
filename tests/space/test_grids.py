#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test grids
"""

import pytest
import random

import stress_fiber as sf

grids = sf.space.grids


grid_list = [grids.HexGrid(n, True) for n in (0, 1, 2, 3)]
grid_list += [grids.HexGrid(n, False) for n in (1, 3)]
coords = ((-2, 3, -1), (0, 0, 0), (1, -1, 0), (-3, 1, 2), (1, 0, -1))

# TODO: Convert to using coords of given size and multiplying them to scale


class TestGrid:
    @pytest.mark.parametrize("grid", grid_list)
    def test__str__(self, grid):
        """Should not return a <blah>"""
        assert not str(grid).startswith("<")

    @pytest.mark.parametrize("grid", grid_list)
    def test_within(self, grid):
        """Origin is always within, 1 off is if size>=1, 10000 isn't"""
        assert grid.within((0, 0, 0))
        if grid.size >= 1:
            assert grid.within((1, -1, 0))
        if grid.size <= 1_000_000:
            assert not grid.within((2_000_000, -1_000_000, -1_000_000))

    @pytest.mark.parametrize("grid", grid_list)
    def test_validate(self, grid):
        """Some locs are more valid than others"""
        assert not grid.validate((0, 1, 0))
        assert not grid.validate((100, 100, 100))
        assert all([grid.validate(loc) for loc in coords])

    @pytest.mark.parametrize("grid", grid_list)
    def test_to_array_indices(self, grid):
        """Created grids should match up to the resolution method"""
        n = grid.size
        for loc in coords:
            if grid.within(loc):
                indices = grid.to_array_indices(loc)
                assert grid.array[indices]["cube"] == loc

    @pytest.mark.parametrize("grid", grid_list)
    def test_neighbors(self, grid):
        """Created grids should match up to the resolution method"""
        loc1 = random.choice(grid.all_entries)["cube"]
        for loc2 in grid.neighbors(loc1):
            assert grid.distance(loc1, loc2) == 1

