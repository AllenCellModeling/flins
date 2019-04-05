# encoding: utf-8
"""
Show me what you got.
CDW 2019
"""

import matplotlib

matplotlib.use("TKAgg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def create_fig():
    fig, ax = plt.subplots(1, 1)
    return fig, ax


def plot_actin(ax, actin):
    """Plot an actin fil onto an axis"""
    raise NotImplementedError


def alpha_actinin(ax, actinin):
    """Plot an α-actinin on the supplied axis"""
    raise NotImplementedError
