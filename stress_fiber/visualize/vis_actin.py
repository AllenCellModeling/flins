# encoding: utf-8
"""
Plot fancy versions of our proteins on matplotlib axes

This was transferred here from the protein classes themselves and is mostly for
storage at this point.
"""

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches


def plot_actin(actin, ax=None, show=False, y=0):
    """Plot an actin filament on an axis"""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.axis("off")
        ax.set(
            xlim=actin.boundaries,
            ylim=(y - 2 * actin._rise, y + 2 * actin._rise),
            aspect=1,
        )
    circ = lambda x, y: matplotlib.patches.CirclePolygon(
        (x, y),
        radius=0.66 * actin._rise,
        resolution=20,
        facecolor="skyblue",
        edgecolor="royalblue",
    )
    # Work through each pair
    for x in [p.x for p in actin.pairs]:
        ax.add_patch(circ(x + actin._rise * 0.1, y + 0.5 * actin._rise))
        ax.add_patch(circ(x - actin._rise * 0.1, y - 0.5 * actin._rise))
    if show:
        plt.show()
    return ax
