#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test hexmath
"""

import stress_fiber.space.hexmath as hm


class TestCube:
    COORDS = ((-2, 3, -1), (0, 0, 0), (1, -1, 0), (-3, 1, 2), (1, 0, -1))

    def test_coordinate_round_trip(self):
        """Make it there and back again"""
        for coord in self.COORDS:
            ax_back = hm.Axial.to_cube(*hm.Cube.to_axial(*coord))
            assert ax_back == coord

    def test_to_array_indices(self):
        """Make sure the indices we are creating in grids match up to the
        resolution method
        """
        for n in (0, 1, 2):
            grid = hm.Cube.create_grid_array(n)
            for coord in self.COORDS:
                if hm.Cube.within_radius(*coord, n):
                    indices = hm.Cube.to_array_indices(*coord, n)
                    print(indices)
                    assert grid[indices]['cube']==coord
