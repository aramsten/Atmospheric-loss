from Calculators.function_solvers import calculate_total_mass_loss
from astropy import units
from astropy.constants import R_sun, L_sun, sigma_sb, M_earth, R_earth
from Data_management import nasa_exoplanet_download
from Data_management.file_manager import save_catalog
import numpy as np
from astropy.table import Table


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
        loss = calculate_total_mass_loss(eta,
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
            colname = f"Loss/pl_mass, t={t:.1f} Gyr"
        elif output == "fraction":
           loss = loss / protoatmosphere_mass_fraction
           colname = f"Loss/{protoatmosphere_mass_fraction}protoatm., t={t:.1f} Gyr"
        else:
            raise ValueError("output must be 'mass' or 'fraction'")
        loss_table[colname] = loss.decompose()
        loss_table[colname].format = ".5f"

    return loss_table

def main():
    catalog = nasa_exoplanet_download.main()

    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr
    R_xuv = 1.2  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01
    description = f"Catalog_mass_loss_for_{end_time[0]}-{end_time[-1]}_Gyr_eta-{eta}_Rxuv-{R_xuv}"
    target = "catalog" # "catalog" or "new"
    output = "mass" # "mass" or "fraction"

    catalog = total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output)
    save_catalog(catalog, "AR", description)

    return catalog

if __name__ == "__main__":
    main()