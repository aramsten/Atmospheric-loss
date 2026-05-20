from Plot_tools.plot_creator import Plot2D_creator, save_plot
from Calculators import function_solvers as fs
from astropy.constants import L_sun, G
from astropy import units as u
import Plot_tools.plot_modifier as pm
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

def plot_cosmic_shoreline_spectral_types(catalog, colname, y_label,names_to_print):
    x_label = "Escape velocity (km/s)"

    x_axis = catalog["v_esc"]
    y_axis = catalog[colname]

    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, spectral_colors = pm.spectral_type_color_coder(catalog, "SpT_PM")

    ax = plot.gca()
    shoreline_position_text = draw_shoreline(ax, catalog, x_axis, y_axis)
    pm.apply_colors(ax, color_per_row, spectral_colors)

    pm.set_point_size_for_names(ax, catalog, "pl_name", names_to_print, color_per_row)
    pm.set_log_axis_base_ten(ax)
    pm.append_text_label(ax, catalog, "pl_name", names_to_print, x_axis, y_axis)

    return plot, shoreline_position_text

def plot_cosmic_shoreline_over_x_with_mass_loss_fraction(catalog, colname, y_label, eta, names_to_print):
    x_label = r"$x = (v_\mathrm{esc}/v_\mathrm{esc\oplus})^3 (\rho / \rho_\oplus)^{0.5}$"

    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
    v_esc = catalog["v_esc"]
    rho = 3*catalog["pl_masse"].to(u.kg)/(4*np.pi*(catalog["pl_rade"].to(u.m))**3)

    v_earth = v_esc.quantity[earth_idx]
    rho_earth = rho[earth_idx]
    x_axis = (v_esc / v_earth) ** 3 * (rho / rho_earth)**0.5
    y_axis = catalog[colname]

    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, spectral_colors = pm.spectral_type_color_coder(catalog, "SpT_PM")

    ax = plot.gca()
    pm.apply_colors(ax, color_per_row, spectral_colors)

    x_vals_raw = x_axis.value if hasattr(x_axis, "value") else x_axis
    y_vals_raw = y_axis.value if hasattr(y_axis, "value") else y_axis
    x_min, x_max = np.nanmin(x_vals_raw), np.nanmax(x_vals_raw)
    y_min, y_max = np.nanmin(y_vals_raw), np.nanmax(y_vals_raw)
    ax.set_xlim(x_min * 0.5, 1e5)
    ax.set_ylim(1e-3, y_max * 2)

    mass_fractions = [2e-4, 2e-3, 1e-2, 2e-1]
    draw_mass_loss_lines(ax, catalog, eta, x_axis, y_axis, mass_fractions)

    pm.set_point_size_for_names(ax, catalog, "pl_name", names_to_print, color_per_row, size_selected=30, size_other=4)
    pm.append_text_label(ax, catalog, "pl_name", names_to_print, x_axis, y_axis, fontweight="bold", auto=True)
    ax.grid(False, which='minor')

    return plot

def plot_cosmic_shoreline_lost_primordial(catalog, colname, y_label, names_to_print):
    x_label = "Escape velocity (km/s)"

    x_axis = catalog["v_esc"]
    y_axis = catalog[colname]

    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, loss_colors = pm.removed_primordal_atmosphere_color_coder(catalog, "Loss/0.01protoatm., star_age")
    ring_colors, spectral_colors = pm.spectral_type_color_coder(catalog, "SpT_PM")

    loss_values = catalog["Loss/0.01protoatm., star_age"]
    color_per_row = [
        ring_color if loss <= 1 else loss_colors["Loss: >0.01 M"]
        for loss, ring_color in zip(loss_values, ring_colors)
    ]

    ax = plot.gca()
    shoreline_position_text = draw_shoreline(ax, catalog, x_axis, y_axis)
    pm.apply_colors(ax, color_per_row, loss_colors)
    pm.set_point_size_for_names(ax, catalog, "pl_name", names_to_print, color_per_row)
    pm.apply_ring_colors(ax, ring_colors, spectral_colors)
    pm.set_log_axis_base_ten(ax)
    pm.append_text_label(ax, catalog, "pl_name", names_to_print, x_axis, y_axis, fontweight="bold", auto=True)

    return plot, shoreline_position_text

def plot_cosmic_shoreline_over_x_with_mass_loss_fraction_lost_primordial(catalog, colname, y_label, eta, names_to_print):
    x_label = r"$x = (v_\mathrm{esc}/v_\mathrm{esc\oplus})^3 (\rho / \rho_\oplus)^{0.5}$"


    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
    v_esc = catalog["v_esc"]
    rho = 3*catalog["pl_masse"].to(u.kg)/(4*np.pi*(catalog["pl_rade"].to(u.m))**3)

    v_earth = v_esc.quantity[earth_idx]
    rho_earth = rho[earth_idx]
    x_axis = (v_esc / v_earth) ** 3 * (rho / rho_earth)**0.5
    y_axis = catalog[colname]

    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, spectral_colors = pm.removed_primordal_atmosphere_color_coder(catalog, "Loss/0.01protoatm., star_age", preserved_color="green")

    ax = plot.gca()
    pm.apply_colors(ax, color_per_row, spectral_colors)

    x_vals_raw = x_axis.value if hasattr(x_axis, "value") else x_axis
    y_vals_raw = y_axis.value if hasattr(y_axis, "value") else y_axis
    x_min, x_max = np.nanmin(x_vals_raw), np.nanmax(x_vals_raw)
    y_min, y_max = np.nanmin(y_vals_raw), np.nanmax(y_vals_raw)
    ax.set_xlim(x_min * 0.5, 1e5)
    ax.set_ylim(1e-3, y_max * 2)

    mass_fractions = [2e-4, 2e-3, 1e-2, 2e-1]
    draw_mass_loss_lines(ax, catalog, eta, x_axis, y_axis, mass_fractions)

    pm.set_point_size_for_names(ax, catalog, "pl_name", names_to_print, color_per_row, size_selected=30, size_other=4)
    pm.append_text_label(ax, catalog, "pl_name", names_to_print, x_axis, y_axis, fontweight="bold", auto=True)
    ax.grid(False, which='minor')

    return plot

def draw_mass_loss_lines(ax, catalog, eta, x_axis, y_axis, fractions):
    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]

    d = catalog["pl_orbsmax"].to(u.m)
    lum = catalog["st_lum"]*L_sun
    lq = catalog["Lxuv/Lbol"]
    t_sat = catalog["t_sat"].to(u.s)
    age = catalog["st_age"].to(u.s)
    gamma = catalog["gamma"]

    d_earth = d[earth_idx]
    lum_sun = lum[earth_idx]
    lq_sun = lq[earth_idx]
    t_sat_sun = t_sat[earth_idx]
    age_sun = age[earth_idx]
    gamma_sun = gamma[earth_idx]

    E_xuv_sun = fs.calculate_total_l_xuv(lum_sun, lq_sun, t_sat_sun, age_sun, gamma_sun)

    v_esc = catalog["v_esc"].to(u.m/u.s)
    rho = 3*catalog["pl_masse"].to(u.kg)/(4*np.pi*(catalog["pl_rade"].to(u.m))**3)
    v_earth = v_esc[earth_idx]
    rho_earth = rho[earth_idx]

    xmin, xmax = ax.get_xlim()
    x_vals = np.logspace(np.log10(xmin), np.log10(xmax), 200)

    ax = plt.gca()

    for f in fractions:
        m = (d_earth ** 2 / (eta * E_xuv_sun)) * (8 * np.pi / (3 * G))**0.5 * f *(v_earth) ** 3 * (rho_earth)**0.5
        y_line = m * x_vals
        ax.plot(x_vals, y_line, linestyle="--", color="black", alpha=0.5, zorder=1)
        x_pos = x_vals[11] if f != 0.0002 else x_vals[15]
        y_pos = y_line[11] if f != 0.0002 else y_line[15]
        ax.text(x_pos, y_pos, rf"$\Delta M = {f}$",
                fontsize=8,
                va="center",
                ha="center",
                rotation =38,
                zorder=2,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1)
                )

def draw_shoreline(ax, catalog, x_axis, y_axis, planetname_column="pl_name", reference="Mars", factor_over_reference=1.1, vmin=1, vmax=100):
    reference_idx = np.where(catalog[planetname_column] == reference)[0][0]
    v_reference = x_axis[reference_idx]
    I_reference = y_axis[reference_idx]

    C = I_reference / (v_reference**4) * factor_over_reference
    y1 = C * vmin**4
    y2 = C * vmax**4

    ax.plot([vmin, vmax], [y1, y2], "--", color="black", linewidth=1.5)
    if factor_over_reference == 1:
        shoreline_position_text = f"through-{reference}"
    else:
        shoreline_position_text = f"at-{factor_over_reference}-{reference}-insolation"

    return shoreline_position_text

def plot_loss(catalog: Table, colname: str, y_label: str, normalize: bool, names_to_print: list) -> plt.Figure:
    """Creates a plot of mass loss or mass loss fraction against escape velocity for a catalog of planets, with the option to normalize the loss to Earth's loss.
     The points are colored by spectral type and labeled with planet names.

    Parameters
    ----------
    catalog : Table
        A catalog containing planet data. Must include columns named "pl_name", "v_esc" and "SpT_PM".
    colname : str
        The name of the column containing the loss data to plot.
    y_label : str
        The label for the y-axis.
    normalize : bool
        Whether to normalize the loss data to the Earth's loss.
    names_to_print : list
        A list of planet names (str) to print labels for.

    Returns
    -------
    plot : matplotlib.pyplot
        The created plot.
    """
    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
    x_axis = catalog["v_esc"]
    if normalize:
        y_axis = catalog[colname] / catalog[colname][earth_idx]

    else:
        y_axis = catalog[colname]

    x_label = "Escape velocity (km/s)"

    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, spectral_colors = pm.spectral_type_color_coder(catalog, "SpT_PM")

    ax = plot.gca()
    pm.apply_colors(ax, color_per_row, spectral_colors)

    xmin, xmax = ax.get_xlim()
    x_vals = np.logspace(np.log10(xmin), np.log10(xmax), 200)
    y_line = [10**0]*len(x_vals)
    ax.plot(x_vals, y_line, linestyle="--",color="black", alpha=0.8, zorder=1)
    x_pos = int(len(x_vals)*0.01)
    y_pos = 10**0
    ax.text(x_pos, y_pos, rf"$10^{{{0}}}$",
            fontsize=8,
            va="center",
            ha="center",
            zorder=2,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1)
            )

    pm.set_point_size_for_names(ax, catalog, "pl_name", names_to_print, color_per_row)
    pm.set_log_axis_base_ten(ax)

    pm.append_text_label(ax, catalog, "pl_name", names_to_print, x_axis, y_axis)

    return plot

def specific_time_plots(catalog, initials, R_xuv, eta, protoatmosphere_mass_fraction, output, loss_plot, normalize_loss, shoreline_plot, t, names_to_print):
    if output == "mass":
        colname_loss = f"Loss/pl_mass, t={t:.1f} Gyr"
    elif output == "fraction":
        colname_loss = f"Loss/{protoatmosphere_mass_fraction}protoatm., t={t:.1f} Gyr"
    else:
        raise ValueError("output must be 'mass' or 'fraction'")
    colname_insol = f"insol_{t}_Gyr"

    if loss_plot:
        if normalize_loss:
            norm_text = "-normalised-to-earth-percent-"
            percent_label = r"\%"
            planet_label = "loss,earth,"
        else:
            norm_text = "-"
            percent_label = ""
            planet_label = "planet"

            proto_mass_frac_text = f"{protoatmosphere_mass_fraction*100 :.0f}%" if output == "fraction" and not normalize_loss else ""

            nominator = rf"$M_{{\mathrm{{loss}}{percent_label}}}$({t}) Gyr"
            denomininator = rf"{proto_mass_frac_text} $M_{{\mathrm{{{planet_label}}}{percent_label}}}$"

            y_label = rf"{nominator} / {denomininator}"
            loss_plot = plot_loss(catalog, colname_loss, y_label, normalize_loss, names_to_print)
            save_plot(loss_plot, initials, f"mass-loss-{proto_mass_frac_text}{norm_text}at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")

    if shoreline_plot:
        y_label = f"Insolation relative to Earth at {t} Gyr"
        cosmic_shoreline_plot, shoreline_position_text = plot_cosmic_shoreline_spectral_types(catalog, colname_insol, y_label, names_to_print)
        save_plot(cosmic_shoreline_plot, initials, f"cosmic_shoreline-spectral_types-{shoreline_position_text}-at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")

def star_age_plots(catalog, initials, R_xuv, eta, protoatmosphere_mass_fraction, output, loss_plot, normalize_loss, shoreline_plot, names_to_print):
    """Plots for mass loss and cosmic shoreline at the stars age"""
    if loss_plot:
        colname_loss = f"Loss/0.01protoatm., star_age"

        if normalize_loss:
            norm_text = "-normalised-to-earth-percent-"
            percent_label = r"\%"
            planet_label = "loss,earth,"
        else:
            norm_text = "-"
            percent_label = ""
            planet_label = "planet"

            proto_mass_frac_text = f"{protoatmosphere_mass_fraction*100 :.0f}%" if output == "fraction" and not normalize_loss else ""

            nominator = rf"$M_{{\mathrm{{loss}}{percent_label}}}$"
            denomininator = rf"{proto_mass_frac_text} $M_{{\mathrm{{{planet_label}}}{percent_label}}}$"

            y_label = rf"{nominator} / {denomininator}"
            loss_plot = plot_loss(catalog, colname_loss, y_label, normalize_loss, names_to_print)
            save_plot(loss_plot, initials, f"mass-loss-{proto_mass_frac_text}{norm_text}at-stars_age-rxuv_factor={R_xuv}-eta={eta}")

    if shoreline_plot:
        y_label = f"Insolation relative to Earth"
        colname_insol = f"insol_star_age"
        cosmic_shoreline_plot, shoreline_position_text = plot_cosmic_shoreline_spectral_types(catalog, colname_insol, y_label, names_to_print)
        save_plot(cosmic_shoreline_plot, initials, f"cosmic_shoreline-spectral_types-{shoreline_position_text}-at-t=star_age-rxuv_factor={R_xuv}-eta={eta}")

def main():
    table_name = "260520_21.27_ST_Catalog_mass_loss_for_0.1-10.0_Gyr_eta-0.1_Rxuv-1.0.ecsv"

    catalog = ascii.read(f"Tables/{table_name}")
    initials = "ST"

    names_to_print = {"Mercury",
                    "Venus",
                    "Earth",
                    "Mars",
                    "Moon", 
                    "Jupiter",
                    "Saturn",
                    "Uranus",
                    "Neptune",
                    "TRAPPIST-1 b",
                    "TRAPPIST-1 c",
                    "TOI-1685 b"}

    R_xuv = catalog.meta["R_xuv"] # dimensionless ratio >= 1
    eta = catalog.meta["eta"] # dimensionless heating efficiency
    protoatmosphere_mass_fraction = catalog.meta["Protoatmosphere mass fraction"]
    output = "fraction" # "mass" or "fraction"
    loss_plot = True
    normalize_loss = False
    shoreline_plot = True
    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr
    end_time = catalog.meta["end_time"]

  #  for t in end_time:
  #      specific_time_plots(catalog, initials, R_xuv, eta, protoatmosphere_mass_fraction, output, loss_plot, normalize_loss, shoreline_plot, t)
    
    star_age_plots(catalog, initials, R_xuv, eta, protoatmosphere_mass_fraction, output, loss_plot, normalize_loss, shoreline_plot,names_to_print)

    plot, shoreline_position_text = plot_cosmic_shoreline_lost_primordial(catalog, "insol_star_age", "Insolation relative to Earth", names_to_print)
    save_plot(plot, initials, f"cosmic_shoreline-lost_primordial-{shoreline_position_text}-at-stars_age-rxuv_factor={R_xuv}-eta={eta}")

    plot = plot_cosmic_shoreline_over_x_with_mass_loss_fraction(catalog, "insol_star_age","Insolation relative to Earth", eta, names_to_print)
    save_plot(plot, initials,f"cosmic_shoreline-at-stars_age-rxuv_factor={R_xuv}-eta={eta}")

    plot = plot_cosmic_shoreline_over_x_with_mass_loss_fraction_lost_primordial(catalog, "insol_star_age","Insolation relative to Earth", eta, names_to_print)
    save_plot(plot, initials, f"cosmic_shoreline-at-stars_age_lost_primordial-rxuv_factor={R_xuv}-eta={eta}")


if __name__ == "__main__":
    main()