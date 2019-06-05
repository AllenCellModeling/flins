
#!/usr/bin/end python
# -*- coding: UTF-8 -*-

"""
Test the tract space and tracts. 
"""

import stress_fiber.space as space

def test_coordinate_round_trip():
    """Make it there and back again"""
    coords = ((-2, 3, -1), (0, 0, 0), (1, -1, 0))
    for coord in coords:
        hex = space.hexmath
        ax_back = hex.Axial.to_cube(*hex.Cube.to_axial(*coord))
        assert ax_back == coord

