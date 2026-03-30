from Calculators.function_solvers import calculate_total_mass_loss
from astropy.constants import R_sun, L_sun, sigma_sb, M_earth, R_earth
from Data_management import nasa_exoplanet_download
from Data_management.file_manager import save_catalog
from Plot_tools.plot_creator import Plot2D_creator, save_plot
import Plot_tools.plot_modifier as pm
import numpy as np
from astropy.table import Table
from astropy import units as u
import matplotlib.pyplot as plt


def total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output, loss_plot, shoreline_plot, normalise_loss):
    """
    target = "new"  --> Makes new table
    target = "catalog" --> Adds columns to catalog

    output = "mass" --> Calculates total mass loss over planet mass
    output = "fraction" --> Calculates total mass loss over 1% of protoatmosphere
    """

    if target == "new":
        loss_table = Table()
        loss_table["pl_name"] = catalog["pl_name"]

    elif target == "catalog":
        loss_table = catalog

    else:
        raise ValueError("target must be 'new' or 'catalog'")

    loss_table.meta["eta"] = eta
    loss_table.meta["R_xuv"] = R_xuv
    loss_table.meta["end_time"] = end_time

    pl_mass = catalog["pl_masse"].to(u.kg)
    pl_rad = catalog["pl_rade"].to(u.m)
    pl_orb = catalog["pl_orbsmax"].to(u.m)
    t_sat = catalog["t_sat"].to(u.s)
    st_lum = catalog["st_lum"] * L_sun

    for t in end_time:
        loss, f_xuv, _ = calculate_total_mass_loss(eta,
                                        pl_mass,
                                        R_xuv * pl_rad,
                                        pl_rad,
                                        st_lum,
                                        catalog["Lxuv/Lbol"],
                                        t_sat,
                                        (t * u.Gyr).to(u.s),
                                        catalog["gamma"],
                                        pl_orb)

        # Normalised mass to planet mass:
        loss = loss / pl_mass

        if output == "mass":
            colname_loss = f"Loss/pl_mass, t={t:.1f} Gyr"
        elif output == "fraction":
           loss = loss / protoatmosphere_mass_fraction
           colname_loss = f"Loss/{protoatmosphere_mass_fraction}protoatm., t={t:.1f} Gyr"
        else:
            raise ValueError("output must be 'mass' or 'fraction'")
        loss_table[colname_loss] = loss.decompose()
        loss_table[colname_loss].format = ".5f"

        colname_insol = f"insol_{t}_Gyr"

        earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
        insolation = f_xuv / f_xuv[earth_idx]

        loss_table[colname_insol] = insolation
        loss_table[colname_insol].format = ".1f"

        plot_loss(loss_table, colname_loss, t, R_xuv, eta, output, protoatmosphere_mass_fraction, normalise_loss) if loss_plot else None
        plot_cosmic_shoreline(loss_table, colname_insol, t, R_xuv, eta) if shoreline_plot else None

    return loss_table

def plot_cosmic_shoreline(catalog, colname, t, R_xuv, eta):
    x_label = "Escape velocity (km/s)"
    y_label = f"Insolation relative to Earth at {t} Gyr"

    x_axis = catalog["v_esc"]
    y_axis = catalog[colname]

    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, spectral_colors = pm.spectral_type_color_coder(catalog, "SpT_PM")

    ax = plot.gca()
    shoreline_position_text = draw_shoreline(ax, catalog, x_axis, y_axis)
    pm.apply_spectral_colors(ax, color_per_row, spectral_colors)

    name_list = {"Mercury", "Venus", "Earth", "Mars"}
    pm.set_point_size_for_names(ax, catalog, "pl_name", name_list, color_per_row)
    pm.set_log_axis_base_ten(ax)
    pm.append_text_label(ax, catalog, "pl_name", name_list, x_axis, y_axis)

    save_plot(plot,"AR",f"cosmic_shoreline-{shoreline_position_text}-at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")

def draw_shoreline(ax, catalog, x_axis, y_axis, planetname_column="pl_name", reference="Mars", factor_over_reference=1.1, vmin=4, vmax=100):
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

def plot_loss(catalog, colname, t, R_xuv, eta, output, proto_mass_frac, normalize):
    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
    x_axis = catalog["v_esc"]
    if normalize:
        y_axis = catalog[colname] / catalog[colname][earth_idx]
        norm_text = "-normalised-to-earth-percent-"
        percent_label = r"\%"
        planet_label = "loss,earth,"
    else:
        y_axis = catalog[colname]
        norm_text = "-"
        percent_label = ""
        planet_label = "planet"

    proto_mass_frac_text = f"{proto_mass_frac}%" if output == "fraction" and not normalize else ""


    nominator = rf"$M_{{\mathrm{{loss}}{percent_label}}}$({t}) Gyr"
    denomininator = rf"{proto_mass_frac_text} $M_{{\mathrm{{{planet_label}}}{percent_label}}}$"

    x_label = "Escape velocity (km/s)"
    y_label = rf"{nominator} / {denomininator}"



    plt.figure()
    plot_creator = Plot2D_creator(x_axis)
    plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    color_per_row, spectral_colors = pm.spectral_type_color_coder(catalog, "SpT_PM")

    ax = plot.gca()
    pm.apply_spectral_colors(ax, color_per_row, spectral_colors)

    name_list = {"Mercury", "Venus", "Earth", "Mars"}
    pm.set_point_size_for_names(ax, catalog, "pl_name", name_list, color_per_row)
    pm.set_log_axis_base_ten(ax)

    pm.append_text_label(ax, catalog, "pl_name", name_list, x_axis, y_axis)

    save_plot(plot, "AR", f"mass-loss-{proto_mass_frac_text}{norm_text}at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")


def main():
    catalog = nasa_exoplanet_download.main()

    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr
    R_xuv = 1.2  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01
    description = f"Catalog_mass_loss_for_{end_time[0]}-{end_time[-1]}_Gyr_eta-{eta}_Rxuv-{R_xuv}"

    target = "catalog" # "catalog" or "new"
    output = "fraction" # "mass" or "fraction"
    loss_plot = True
    normalise_loss = False
    shoreline_plot = True
    save_to_catalog = False

    catalog = total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output, loss_plot, shoreline_plot, normalise_loss)
    save_catalog(catalog, "AR", description) if save_to_catalog else None
    return catalog

if __name__ == "__main__":
    main()