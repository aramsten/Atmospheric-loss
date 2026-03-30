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


def spectral_type_color_coder(catalog, colname="SpT_PM"):
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

def apply_spectral_colors(ax, color_per_row, spectral_colors):
    points = ax.collections[0]
    points.set_facecolors(color_per_row)

    # Legend-proxy-punkter
    for spt, col in spectral_colors.items():
        ax.scatter([], [], color=col, s=30, label=spt)

    ax.legend(loc="best")

def set_point_size_for_names(ax, catalog, name_column, names, color_per_row,
                             size_selected=20, size_other=5,
                             alpha_selected=0.9, alpha_other=0.7,
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
                      color="green", fontsize=8):
    for i, name in enumerate(catalog[colname]):
        if name not in names:
            continue
        x = x_axis[i]
        y = y_axis[i]
        ax.text(x*x_offset,y*y_offset,name,fontsize=fontsize,color=color,ha="left",va="bottom")
