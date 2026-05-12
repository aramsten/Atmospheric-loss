from matplotlib.table import Table
import numpy as np
from matplotlib.ticker import ScalarFormatter, LogLocator

def set_log_axis_base_ten(ax, axis_to_modify="x"):
    if axis_to_modify == "x":
        axis = ax.xaxis
    elif axis_to_modify == "y":
        axis = ax.yaxis
    else:
        raise ValueError("axis_to_modify must be 'x' or 'y'")

    axis.set_major_formatter(ScalarFormatter())
    axis.set_major_locator(LogLocator(base=10))


def spectral_type_color_coder(catalog: Table, colname="SpT_PM") -> list[str] | dict[str, str]:
    """Assign colors to spectral types in the catalog based on the specified column.

    Parameters
    ----------
    catalog : Table
        The catalog containing the spectral type information.
    colname : str, optional
        The name of the column containing the spectral types (default is "SpT_PM").

    Returns
    -------
    colors : list[str]
        A list of colors corresponding to each row in the catalog based on spectral type.
    spectral_colors : dict[str, str]
        A dictionary mapping spectral types to their corresponding colors.
    """

    colors = []
    spectral_colors = {
        "F": "blue",
        "G": "green",
        "K": "orange",
        "M": "red"}

    for spt in catalog[colname]:
        letter = str(spt)[0]
        color = spectral_colors.get(letter, "gray")
        colors.append(color)

    return colors, spectral_colors

def removed_primordal_atmosphere_color_coder(
        catalog: Table,
        colname="Loss/0.01protoatm., star_age",
        removed_color = "lightgrey",
        preserved_color = "black") \
        -> list[str] | dict[str, str]:
    """Assign colors to planets based on whether they have lost their primordial atmosphere or not.
    
    Parameters
    ----------
    catalog : Table
        The catalog containing the atmospheric loss information.
    colname : str, optional
        The name of the column containing the atmospheric loss values (default is "Loss/0.01protoatm., star_age").

    Returns
    -------
    colors : list[str]
        A list of colors corresponding to each row in the catalog based on atmospheric loss.
    loss_colors : dict[str, str]
        A dictionary mapping atmospheric loss states to their corresponding colors.
    """
    colors = []
    loss_colors = {
        "Loss: >0.01 M": removed_color,
        "Loss: <0.01 M": preserved_color
    }

    for loss in catalog[colname]:
        if loss > 1:  # Assuming positive values indicate removed atmosphere
            colors.append(loss_colors["Loss: >0.01 M"])
        else:
            colors.append(loss_colors["Loss: <0.01 M"])
    return colors, loss_colors

def apply_colors(ax, color_per_row, colors_dir):
    """Apply colors to the points in the plot and create a legend.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to which the colors will be applied.
    color_per_row : list[str]
        A list of colors for each row in the catalog.
    colors_dir : dict[str, str]
        A dictionary mapping types to their corresponding colors.
    """

    points = ax.collections[0]
    points.set_facecolors(color_per_row)

    # Legend-proxy-punkter
    for spt, col in colors_dir.items():
        ax.scatter([], [], color=col, s=30, label=spt)

    ax.legend(loc="best")

def apply_ring_colors(ax, color_per_row, colors_dir):
    """Apply colors to the edges of the points in the plot.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to which the colors will be applied.
    color_per_row : list[str]
        A list of colors for each row in the catalog.
    colors_dir : dict[str, str]
        A dictionary mapping types to their corresponding colors.
    """
    points = ax.collections[0]
    points.set_edgecolors(color_per_row)

        # Legend-proxy-punkter
    for spt, col in colors_dir.items():
        ax.scatter([], [], color=col, s=30, label=spt)

    ax.legend(loc="best")

def set_point_size_for_names(ax, catalog, name_column, names, color_per_row,
                             size_selected=20, size_other=5,
                             alpha_selected=1, alpha_other=0.6,
                             edge_color_selected="black",
                             edge_width_selected=0.4):
    points = ax.collections[0]

    is_selected = np.array([name in names for name in catalog[name_column]])

    n = len(catalog)
    sizes = np.zeros(n)
    alpha = np.zeros(n)
    edgecols = []

    for i in range(n):
        if is_selected[i]:
            sizes[i] = size_selected
            edgecols.append(edge_color_selected)
            alpha[i] = alpha_selected
        else:
            sizes[i] = size_other
            edgecols.append(color_per_row[i])
            alpha[i] = alpha_other

    points.set_sizes(sizes)
    points.set_edgecolors(edgecols)
    points.set_alpha(alpha)
    points.set_linewidth(edge_width_selected)

def append_text_label(ax, catalog, colname, names, x_axis, y_axis,
                      x_offset=1.05, y_offset=1.02,
                      color="black", fontweight="normal", fontsize=8, auto=False):
    texts = []

    if auto:
        x_offset = 1
        y_offset = 1
    for i, name in enumerate(catalog[colname]):

        if name not in names:
            continue
        x = x_axis[i]
        y = y_axis[i]

        # Nedan omvandlar x,y till float omm x,y är quantity
        if hasattr(x, "to_value"):
            x = x.to_value()
        if hasattr(y, "to_value"):
            y = y.to_value()

        texts.append(ax.text(
                x*x_offset,
                y*y_offset,
                name,
                va="center",
                ha="left",
                fontsize=fontsize,
                color=color,
                fontweight=fontweight,
                ))

    if auto:
        adjust_text(
            texts,
            x=x_axis,
            y=y_axis,
            ax=ax,

            # These five parameters should change the position of the text/length of arrow
            # It seems to be working rather poorly for our log-plots
            # Change with caution
            force_points=(0.5, 1.0),
            force_text=(0.5, 1.0),
            expand_points=(0.5, 1.0),
            axpand_text=(0.5, 1.0),
            min_arrow_len=15,

            arrowprops=dict(
                arrowstyle='->',
                color='grey',
                lw=0.8,
                # Shrink should change the starting point of the arrow
                # Works rather poorly for our log-plot
                shrinkA=8,
                shrinkB=1
            )
        )