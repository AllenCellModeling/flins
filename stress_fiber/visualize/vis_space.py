# encoding: utf-8
"""
Plot things about our hex spaces and their tracts on matplotlib axes
"""

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches
from ..space.hexmath import Cube


def plot_tractspace(ts, callback=None, show=False):
    """Plot the tracts and what's in them, optionally using a callback"""
    # Set up callback
    if callback is None:
        callback = lambda t: t.loc  # default to cube coordinates
    # Set up figure
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.axis("off")
    limit = 2 * ts.size + 1
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect=1)
    # Work through each tract
    for tract in ts.all_tracts:
        x, y = Cube.to_cart(*tract.loc)
        hex = matplotlib.patches.RegularPolygon(
            (x, y), numVertices=6, radius=1, facecolor="White", edgecolor="k"
        )
        ax.add_patch(hex)
        ax.text(x, y, str(callback(tract)), ha="center", va="center")
    if show:
        plt.show()
    return ax
