# encoding: utf-8
"""
Support for rendering the tracts (and proteins within them) at the current
timestep as an SVG. When changing this, [a list of named
colors](https://www.december.com/html/spec/colorsvg.html) is useful.
"""

import stress_fiber as sf

import numpy as np
import svgwrite


# Can globally set the style of the render via css
CSS_STYLES = """
    .background { background-color: oldlace; }
    .actin { stroke: cornflowerblue; stroke-width: 3px; }
    .anchor { stroke: tan; fill: tan;}
    .actinin { stroke: limegreen; fill: limegreen; stroke-width: 2px; }
    .tract { stroke: darkslategray; fill: whitesmoke; stroke-width: 3px; }
    .fade {opacity: 0.3;}
"""


def _plot_anchor(anchor, group, y, params):
    """Plot anc as part of group at y"""
    x = anchor.x * params["xm"]
    group.add(svgwrite.shapes.Circle((x, y), 4, class_="anchor"))
    return


def _plot_actinin_head(actinin_head, group, y, params):
    """Plot actinin as part of group from y[0] to y[1]"""
    x_i = actinin_head.x * params["xm"]
    if actinin_head.other_head.bs.bound:
        x_f = actinin_head.other_head.x * params["xm"]
    else:
        x_f = x_i
    group.add(svgwrite.shapes.Line((x_i, y[0]), (x_f, y[1]), class_="actinin"))
    group.add(svgwrite.shapes.Circle((x_i, y[0]), 2, class_="actinin"))
    return


def _plot_actinin(actinin, group, y, params):
    """Plot an actinin as part of group within y_lim"""
    y *= params["ym"]
    x = actinin.x * params["xm"]
    group.add(svgwrite.shapes.Ellipse((x, y), (3, 1)))


def _plot_actin(actin, group, y, params):
    """Plot actin and affiliated proteins as part of group within y_lim"""
    # Locally load params
    xm = params["xm"]
    ym = params["ym"]
    y_span = params["y_span"]
    # Find limits and convert to user units
    y_edge = np.round(y / y_span) * y_span * ym  # end links at top or bottom
    x = np.multiply(xm, actin.boundaries)
    y *= ym
    # Plot actin
    group.add(svgwrite.shapes.Line((x[0], y), (x[1], y)))
    # Plot associated mols
    for pair in actin.pairs:
        if not pair.bs.bound:
            continue
        mol = pair.bs.linked
        if type(mol) is sf.alpha_actinin.ActininHead:
            _plot_actinin_head(mol, group, (y, y_edge), params)
        elif type(mol) is sf.anchor.Anchor:
            _plot_anchor(mol, group, y, params)
    return


def _plot_tract(dwg, tract, tract_i, params):
    """Plot tract in dwg offset to loc for tract_i"""
    # Locally load params
    ym = params["ym"]
    y_span = params["y_span"]
    y_sep = params["y_sep"]
    # Find tract y loc
    tract_y_offset = ym * (tract_i * y_span + (tract_i + 0.5) * y_sep)
    # Create tract group and plot outline
    group = dwg.add(dwg.g())
    group.add(
        svgwrite.shapes.Rect(size=(dwg.attribs["width"], ym * y_span), class_="tract")
    )
    group.translate(0, tract_y_offset)
    return group


def plot_world(world, params={}):
    """Plot a world as a flat, unrolled, svg
    Paramters
    ---------
    world: stress_fiber.construct.World
        World to render as svg
    """
    ## Visualization parameters
    # Manage parameters for passing to sub-plotting
    defaults = {
        "y_span": 50,  # how tall we treat tracts as being
        "y_sep": 4,  # distance between tracts
        "xm": 2,  # unit to pixel multiplier
        "ym": 2,  # unit to pixel multiplier
    }
    for key, val in defaults.items():
        if not key in params:
            params[key] = val
    # Locally used params
    y_span = params["y_span"]
    y_sep = params["y_sep"]
    xm = params["xm"]
    ym = params["ym"]
    n_tracts = len(world.tractspace.all_tracts)  # n to plot
    x_span = world.tractspace.span  # length of tracts
    y_tot = n_tracts * (y_sep + y_span)  # calculated total height
    ## Plotting
    # Create drawing
    dwg = svgwrite.Drawing(size=(xm * x_span, ym * y_tot), class_="background")
    dwg.defs.add(dwg.style(CSS_STYLES))  # lets us use class defs from CSS_STYLES
    # Create each tract and plot contents
    for i, tract in enumerate(world.tractspace.all_tracts):
        tract_group = _plot_tract(dwg, tract, i, params)
        # Plot unbound actinin
        n_actinin = len(tract.mols["actinin"])
        actinin_group = tract_group.add(svgwrite.container.Group(class_="actinin fade"))
        for j, actinin in enumerate(tract.mols["actinin"]):
            y_actinin = y_span * (j + 1) / (n_actinin + 1)
            if not actinin.bound:
                _plot_actinin(actinin, actinin_group, y_actinin, params)
        # Plot Actin and associated
        n_actin = len(tract.mols["actin"])
        actin_group = tract_group.add(svgwrite.container.Group(class_="actin"))
        for k, act in enumerate(tract.mols["actin"]):
            y_act = y_span * (k + 1) / (n_actin + 1)
            _plot_actin(act, actin_group, y_act, params)
    return dwg
