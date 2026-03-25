from Calculators.function_solvers import calculate_total_mass_loss
from astropy import units
from astropy.constants import R_sun, L_sun, sigma_sb, M_earth, R_earth
from Data_management import nasa_exoplanet_download
from Data_management.file_manager import save_catalog
from Plot_tools.plot_creator import Plot2D_creator, save_plot
import numpy as np
from astropy.table import Table
from astropy import units as u


def total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output):
    """
    target = "new"  --> Makes new table
    target = "catalog" --> Adds columns to catalog
    output = "mass" --> Calculates total mass loss normalised to M_sun
    output = "fraction" --> Calculates total mass loss normalised to 1% of protoatmosphere
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

    for t in end_time:
        loss, f_xuv, _ = calculate_total_mass_loss(eta,
                                        (catalog["pl_masse"]).to(units.kg),
                                        (R_xuv * catalog["pl_rade"]).to(units.m),
                                        (catalog["pl_rade"]).to(units.m),
                                        catalog["st_lum"]*L_sun,
                                        catalog["Lxuv/Lbol"],
                                        catalog["t_sat"].to(units.s),
                                        (t * units.Gyr).to(units.s),
                                        catalog["gamma"],
                                        catalog["pl_orbsmax"].to(units.m))

        # Normalised mass to planet mass:
        loss = loss / (catalog["pl_masse"]).to(units.kg)
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
        insolation = f_xuv/f_xuv[0]
        loss_table[colname_insol] = insolation
        loss_table[colname_insol].format = ".1f"

        name_loss = "mass-loss-1percent-proto-norm-to-earth"
        x_label_loss = "Escape velocity (km/s)"
        y_label_loss = rf"$M_{{\mathrm{{loss}}}} \, / \, M_{{\mathrm{{loss}},\oplus}}$, % at {t} Gyr"
        #plot_loss(loss_table, colname_loss, t, R_xuv, eta, x_label_loss,y_label_loss, name_loss)

        name_f_xuv = "insolation_rel_to_earth"
        x_label_f_xuv = "Escape velocity (km/s)"
        y_label_f_xuv = f"Insolation relative to Earth, % at {t} Gyr"
        plot_loss(loss_table, colname_insol, t, R_xuv, eta, x_label_f_xuv, y_label_f_xuv, name_f_xuv)

    return loss_table

def plot_loss(catalog, colname, t, R_xuv, eta, x_label, y_label, name):
    plot_creator = Plot2D_creator(catalog["v_esc"])
    norm_loss_to_earth_loss = catalog[colname]/catalog[colname][0]


    plot_creator.append_y_axis(norm_loss_to_earth_loss)
    plot = plot_creator.create_2D_plot(plot_scatter=True,x_label=x_label, y_label=y_label, label=catalog["pl_name"],x_logscale=True, y_logscale=True, view_legend=False)

    from matplotlib.ticker import ScalarFormatter, LogLocator
    ax = plot.gca()
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_major_locator(LogLocator(base=10))
    earth_x = catalog["v_esc"][0]
    earth_y = norm_loss_to_earth_loss[0]
    ax.scatter(earth_x, earth_y, s=80, color="blue", edgecolor="white", zorder=5)
    ax.text(earth_x * 1.05, earth_y, "Earth", color="blue", fontsize=12)

    save_plot(plot,
              "AR",
              f"{name}-at-t={t}Gyr-rxuv_factor={R_xuv}-eta={eta}")


def main():
    catalog = nasa_exoplanet_download.main()

    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr
    R_xuv = 1.2  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01
    description = f"Catalog_mass_loss_for_{end_time[0]}-{end_time[-1]}_Gyr_eta-{eta}_Rxuv-{R_xuv}"
    target = "catalog" # "catalog" or "new"
    output = "fraction" # "mass" or "fraction"

    catalog = total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output)
    #save_catalog(catalog, "AR", description)



    return catalog

if __name__ == "__main__":
    main()