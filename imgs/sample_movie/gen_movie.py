# encoding: utf-8

"""
Create a movie of a simple model.
"""

import stress_fiber as sf
import stress_fiber.visualize

import os
import subprocess
import shutil

## Supporting functions
to_svg = lambda dwg, f: sf.visualize.svg.save(dwg, f+".svg")
def to_png(outputfn, zoom=1.0):
    call = "rsvg-convert -o %s.png -z %0.2f %s.svg"%(outputfn, zoom, outputfn)
    subprocess.call(call, shell=True)
    return
def to_mp4(framerate, fnformat="%03d", outputfn="out.mp4"):
    call = "ffmpeg -y -r %f -i %s.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p %s"%(
        framerate, fnformat, outputfn)
    subprocess.call(call, shell=True)
    return

## Create world
world = sf.construct.create_test_world(1, 2000, 10, 100)

## Create space to save world
tmp_dir = "./tmp/"
if not os.path.isdir(tmp_dir):
    os.mkdir(tmp_dir)

## Run world and save SVGs and PNGs
for i in range(100):
    print(i)
    dwg = sf.visualize.flat_render.plot_world(world)
    fn = os.path.join(tmp_dir, "%03i"%i)
    to_svg(dwg, fn)
    to_png(fn, 1)
    world.step()

## Convert to movie and trash tmp
to_mp4(2, tmp_dir+"%03d")
shutil.rmtree(tmp_dir)
