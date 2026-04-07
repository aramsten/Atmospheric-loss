from Calculators.function_solvers import calculate_total_mass_loss
from astropy.constants import R_sun, L_sun, sigma_sb, M_earth, R_earth
from Data_management import nasa_exoplanet_download
from Data_management.file_manager import save_catalog
import Plot_tools.plot_modifier as pm
import numpy as np
from astropy.table import Table
from astropy import units as u


def total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output):
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

    return loss_table

def main():
    catalog = nasa_exoplanet_download.main()
    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr

    output = "fraction" # "mass" or "fraction"
    R_xuv = 1.2  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01
    
    description = f"Catalog_mass_loss_for_{end_time[0]}-{end_time[-1]}_Gyr_eta-{eta}_Rxuv-{R_xuv}"


    target = "catalog" # "catalog" or "new"
    catalog = total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output)
    save_catalog(catalog, "ST", description)
    
if __name__ == "__main__":
    main()