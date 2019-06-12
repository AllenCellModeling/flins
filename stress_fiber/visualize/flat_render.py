# encoding: utf-8
"""
Support for rendering the tracts (and proteins within them) at the current
timestep as an SVG. When changing this, [a list of named
colors](https://www.december.com/html/spec/colorsvg.html) is useful.
"""

import stress_fiber as sf

import numpy as np
import svgwrite


# Can globally set the style of the render via CSS
CSS_STYLES = """
    .background { background-color: oldlace; }
    .actin { stroke: cornflowerblue; stroke-width: 3px; }
    .anchor { stroke: tan; fill: tan;}
    .actinin { stroke: limegreen; fill: limegreen; stroke-width: 2px; }
    .tract { stroke: darkslategray; fill: whitesmoke; stroke-width: 3px; }
    .clear { fill-opacity: 0.0; }
    .fade { opacity: 0.3; }
"""


def _plot_anchor(anchor, params):
    """Plot anc as part of group at y"""
    x = anchor.x * params["xm"]
    y = anchor.bs.linked.filament.__y * params["ym"]
    group = anchor.bs.linked.filament.tract.__groups['anchor']
    group.add(svgwrite.shapes.Circle((x, y), 4, class_="anchor"))
    return


def _plot_actinin(actinin, y, params):
    """Plot an actinin as part of group within y_lim"""
    if not actinin.bound:  # none bound
        y *= params["ym"]
        x = actinin.x * params["xm"]
        group = actinin.tract.__groups["unbound_actinin"]
        group.add(svgwrite.shapes.Ellipse((x, y), (3, 1)))
        return
    if sum([h.bs.bound for h in actinin.heads]) == 2:  # both bound
        _plot_actinin_body(actinin, params)
    for head in actinin.heads:
        _plot_actinin_head(head, params)
    return


def _plot_actinin_head(head, params):
    """Plot actinin head on actin"""
    if not head.bs.bound:
        return
    x = head.x * params["xm"]
    y = head.bs.linked.filament.__y * params["ym"]
    group = head.bs.linked.filament.tract.__groups["bound_actinin"]
    group.add(svgwrite.shapes.Circle((x, y), 2, class_="actinin"))
    return


def _plot_actinin_body(actinin, params):
    """Only plot body in case where both heads are bound"""
    tracts = [head.bs.linked.filament.tract for head in actinin.heads]
    if tracts[0] == tracts[1]:  # same tract
        x = [head.x * params["xm"] for head in actinin.heads]
        y = [h.bs.linked.filament.__y * params["ym"] for h in actinin.heads]
        group = tracts[0].__groups["bound_actinin"]
        group.add(svgwrite.shapes.Line((x[0], y[0]), (x[1], y[1]), class_="actinin"))
    else:  # different tracts
        for head, tract in zip(actinin.heads, tracts):
            x_i = head.x * params["xm"]
            x_f = head.other_head.x * params["xm"]
            ym, y_span = params["ym"], params["y_span"]
            y = head.bs.linked.filament.__y 
            y_i = y * ym
            y_f = np.round(y / y_span) * y_span * ym  # end links at top or bottom
            group = tract.__groups["bound_actinin"]
            group.add(svgwrite.shapes.Line((x_i, y_i), (x_f, y_f), class_="actinin"))
        return 


def _plot_actin(actin, params):
    """Plot actin and affiliated proteins as part of group within y_lim"""
    # Locally load params
    xm = params["xm"]
    ym = params["ym"]
    y_span = params["y_span"]
    y = actin.__y
    group = actin.tract.__groups["actin"]
    # Find limits and convert to user units
    y_edge = np.round(y / y_span) * y_span * ym  # end links at top or bottom
    x = np.multiply(xm, actin.boundaries)
    y *= ym
    # Plot actin
    group.add(svgwrite.shapes.Line((x[0], y), (x[1], y)))
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
    tract.__groups = {
        "self": group,
        "actin": group.add(svgwrite.container.Group(class_="actin")),
        "anchor": group.add(svgwrite.container.Group(class_="anchor")),
        "bound_actinin": group.add(svgwrite.container.Group(class_="actinin")),
        "unbound_actinin": group.add(svgwrite.container.Group(class_="actinin fade")),
    }
    group.add(
        svgwrite.shapes.Rect(size=(dwg.attribs["width"], ym * y_span),
                             class_="clear tract")
    )
    return group


def plot_world(world, params={}):
    """Plot a world as a flat, unrolled, svg
    Parameters
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
    # Tracts
    for i, tract in enumerate(world.tractspace.all_tracts):
        _plot_tract(dwg, tract, i, params)
    # Actins
    for tract in world.tractspace.all_tracts:
        actins = tract.mols["actin"]
        n_actin = len(actins)
        for k, actin in enumerate(actins):
            actin.__y = y_span * (k + 1) / (n_actin + 1)
            _plot_actin(actin, params)
    # Alpha-actinins
    for tract in world.tractspace.all_tracts:
        actinins = tract.mols["actinin"]
        n_actinin = len(actinins)
        for j, actinin in enumerate(actinins):
            y_actinin = y_span * (j + 1) / (n_actinin + 1)
            _plot_actinin(actinin, y_actinin, params)
    # Anchors
    for tract in world.tractspace.all_tracts:
        if "anchor" in tract.mols:
            for anchor in tract.mols["anchor"]:
                _plot_anchor(anchor, params)
    return dwg
