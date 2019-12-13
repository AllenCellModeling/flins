#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Test movie maker
"""

import pytest

import os
import tempfile
import flins as fl
import flins.visualize


dirs = [None, tempfile.mkdtemp(), "non-existing_test_movie_dir"]


@pytest.mark.parametrize("dir", dirs)
def test_MovieGen(dir):
    world = fl.construct.create_test_world(0, 2000, 1, 2, 2)
    quiet = False if dir is None else True
    moviegen = fl.visualize.movie.MovieGen(temp_dir=dir, quiet=quiet)
    moviegen.add_world(world)
    world.step()
    moviegen.add_world(world)
    # Test existence of pngs and svgs
    svgpath = os.path.join(moviegen.temp_dir, "000001.svg")
    pngpath = os.path.join(moviegen.temp_dir, "000001.png")
    assert os.path.exists(svgpath), "first svg doesn't exist"
    assert os.path.exists(svgpath), "first svg doesn't exist"
    # Write movie and test existance
    moviegen.write_movie()
    assert os.path.exists(moviegen.outname), "movie doesn't exist"
    # Delete and test removals
    os.remove("out.mp4")
    moviegen.clean_up()
    assert not os.path.exists(svgpath), "first svg not deleted"
    assert not os.path.exists(pngpath), "first png not deleted"
