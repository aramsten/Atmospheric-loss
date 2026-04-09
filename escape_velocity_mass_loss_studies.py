from Plot_tools.plot_creator import Plot2D_creator, save_plot
import Plot_tools.plot_modifier as pm
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

def plot_cosmic_shoreline(catalog, colname, y_label):
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
    pm.apply_spectral_colors(ax, color_per_row, spectral_colors)

    name_list = {"Mercury", "Venus", "Earth", "Mars"}
    pm.set_point_size_for_names(ax, catalog, "pl_name", name_list, color_per_row)
    pm.set_log_axis_base_ten(ax)
    pm.append_text_label(ax, catalog, "pl_name", name_list, x_axis, y_axis)

    return plot, shoreline_position_text

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

    return shoreline_position_text, shoreline_position_text

def plot_loss(catalog, colname, y_label, normalize):
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
    pm.apply_spectral_colors(ax, color_per_row, spectral_colors)

    name_list = {"Mercury", "Venus", "Earth", "Mars"}
    pm.set_point_size_for_names(ax, catalog, "pl_name", name_list, color_per_row)
    pm.set_log_axis_base_ten(ax)

    pm.append_text_label(ax, catalog, "pl_name", name_list, x_axis, y_axis)

    return plot


def main():
    table_name = "260409_22.25_ST_Catalog_mass_loss_for_0.1-10.0_Gyr_eta-0.1_Rxuv-1.2.csv"
    catalog = ascii.read(f"Tables/{table_name}")
    initials = "ST"

    R_xuv = 1.2  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01
    output = "fraction" # "mass" or "fraction"
    loss_plot = True
    normalize_loss = False
    shoreline_plot = True
    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr


    for t in end_time:
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
                loss_plot = plot_loss(catalog, colname_loss, y_label, normalize_loss)
                save_plot(loss_plot, initials, f"mass-loss-{proto_mass_frac_text}{norm_text}at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")

        if shoreline_plot:
            y_label = f"Insolation relative to Earth at {t} Gyr"
            cosmic_shoreline_plot, shoreline_position_text = plot_cosmic_shoreline(catalog, colname_insol, y_label)
            save_plot(cosmic_shoreline_plot, initials, f"cosmic_shoreline-{shoreline_position_text}-at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")
    
    if loss_plot:
        colname_loss = f"Loss/0.01protoatm., star_lifetime"

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
            loss_plot = plot_loss(catalog, colname_loss, y_label, normalize_loss)
            save_plot(loss_plot, initials, f"mass-loss-{proto_mass_frac_text}{norm_text}at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")

    if shoreline_plot:
        y_label = f"Insolation relative to Earth during the planet's star lifetime"
        colname_insol = f"insol_star_lifetime"
        cosmic_shoreline_plot, shoreline_position_text = plot_cosmic_shoreline(catalog, colname_insol, y_label)

if __name__ == "__main__":
    main()