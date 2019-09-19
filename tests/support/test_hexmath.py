#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test hexmath
"""

import stress_fiber.support.hexmath as hm


class TestCube:
    COORDS = ((-2, 3, -1), (0, 0, 0), (1, -1, 0), (-3, 1, 2), (1, 0, -1))

    def test_coordinate_round_trip(self):
        """Make it there and back again"""
        for coord in self.COORDS:
            ax_back = hm.axial.to_cube(*hm.cube.to_axial(*coord))
            assert ax_back == coord
