# encoding: utf-8
"""
Support for movie generation
"""

import os
import subprocess
import shutil
import tempfile
import tqdm
import IPython.display
from . import svg
from . import flat_render


class MovieGen:
    """Generate a movie from a set of steps in the world or SVGs"""

    def __init__(
        self, outname="out", temp_dir=None, fps=25, zoom=(1.0, 1.0), quiet=False
    ):
        """Support the creation of movies from runs

        Parameters
        ---------
        outname: str (out)
            Relative path for the output file
        temp_dir: str
            Dir where we'll do our work, deleted on exit if not specified
        fps: int
            How many frames per sec for resulting movie
        zoom: tuple of floats
            x and y zoom levels
        quiet: boolean
            don't show current status if True
        """
        # Create the dir where we'll do the processing
        if temp_dir is None:
            self.temp_dir = tempfile.mkdtemp()
            self.temp_delete = True
        else:
            self.temp_dir = temp_dir
            if not os.path.exists(temp_dir):
                os.makedirs(temp_dir)
                self.temp_delete = True
            else:
                self.temp_delete = False
        # Create the lists to store our processes
        self.svgs = []
        self.pngs = []
        self.outname = outname if outname.endswith(".mp4") else outname + ".mp4"
        # Save params
        self._fps = fps
        self._zoom = zoom
        self._quiet = quiet
        self._save_renders = False  # secret way to save files
        self._movie_written = False

    def clean_up(self):
        """Delete the temporary renders and the temp dirs if they were
        automatically created
        """
        if self.temp_delete:
            shutil.rmtree(self.temp_dir, True)
        elif not self._save_renders:
            for fn in self.svgs + self.pngs:
                try:
                    os.remove(os.path.join(self.temp_dir, fn))
                except FileNotFoundError:
                    pass  # we don't care if the files were already cleaned

    def __del__(self):
        """Call clean up if not already done"""
        self.clean_up()

    def add_world(self, world):
        """Add a world and render it as an svg"""
        dwg = flat_render.plot_world(world)
        self.add_svg(dwg)

    def add_svg(self, dwg):
        """Render an svg to the temp dir"""
        assert len(self.svgs) < 1_000_000
        filename = os.path.join(self.temp_dir, "%06i.svg" % (len(self.svgs) + 1))
        svg.save(dwg, filename)
        self.svgs.append(filename)

    def write_movie(self):
        """Convert svgs to pngs and then to an mp4 via ffmpeg"""
        # Writing svgs to pngs
        self.pngs = []
        if self._quiet:
            svg_iter = self.svgs
        else:
            svg_iter = tqdm.tqdm(self.svgs, "svg->png", leave=False)
        for svg_fn in svg_iter:
            fn = os.path.splitext(svg_fn)[0]
            xz, yz = self._zoom
            call = "rsvg-convert -o %s.png -x %.4f -y %.4f %s.svg" % (fn, xz, yz, fn)
            subprocess.run(call, shell=True)
            self.pngs.append(fn + ".png")
        # Convert pngs to mp4
        call = "ffmpeg -y -r %f " % self._fps
        call += "-i %s " % os.path.join(self.temp_dir, "%06d.png")
        call += "-c:v libx264 -vf fps=25 -pix_fmt yuv420p %s " % self.outname
        if self._quiet:
            call += "-loglevel panic"
        subprocess.run(call, shell=True)
        self._movie_written = True

    def show(self):
        """Use ipython to display movie, intended for use in notebooks"""
        if not self._movie_written:
            self.write_movie()
        IPython.display.Video(self.outname)
