# encoding: utf-8
"""
Support for saving and displaying SVGs
"""

import os
import subprocess
import uuid
from IPython.display import SVG


def save(dwg, filename):
    """Save an SVGA to a local file

    Parameters
    ----------
    dwg: SVG write.Drawing
        SVG to render as a string and save to disk
    filename: str
        Local filename to save to. Will append '.svg' if not present.
    """
    if not filename.endswith("svg"):
        filename += ".svg"
    with open(filename, "w") as fileout:
        fileout.write(dwg.tostring())


def display(dwg):
    """Save to a temporary file and open (only on MacOS)

    Parameters
    ----------
    dwg: SVG write.Drawing
        SVG to render to disk and open in Safari
    """
    filename = os.environ["TMPDIR"] + str(uuid.uuid1()) + ".svg"
    with open(filename, "w") as tmpfile:
        tmpfile.write(dwg.tostring())
    subprocess.run(["open", "-aSafari", filename])


def jupyter(dwg):
    """Show the drawing in jupyter

    Parameters
    ----------
    dwg: SVG write.Drawing
        SVG to render via IPython.display
    """
    return SVG(dwg.tostring())
