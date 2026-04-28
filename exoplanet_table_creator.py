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

    # Saves meta-data to for calculated parameters catalog for future referencing
    loss_table.meta["eta"] = eta
    loss_table.meta["R_xuv"] = R_xuv
    loss_table.meta["end_time"] = end_time
    catalog.meta["Protoatmosphere mass fraction"] = protoatmosphere_mass_fraction


    for t in end_time:
        loss_table = add_loss_columns_specific_time(loss_table,catalog, t, R_xuv, eta, protoatmosphere_mass_fraction, output)
    
    loss_table = add_loss_colums_star_age(loss_table, catalog, R_xuv, eta, protoatmosphere_mass_fraction, output)

    return loss_table

def add_loss_columns_specific_time(loss_table: Table, catalog: Table, end_time: float, R_xuv: float, eta: float, protoatmosphere_mass_fraction: float, output: str) -> Table:
    """Add a colum to the table with the mass loss for a specific time

    Parameters
    ----------
    loss_table: Table
        the table to which the column will be added
    catalog: Table
        the original catalog with the planet data
    end_time: float
        the time at which to calculate the mass loss
    R_xuv: float
        the ratio of the XUV radius to the planet radius
    eta: float
        the heating efficiency
    protoatmosphere_mass_fraction: float
        the fraction of the planet's mass that is in the protoatmosphere
    output: str
        the type of output to calculate ("mass" or "fraction")

    Returns
    -------
    loss_table: Table
        the updated table with the added column"""
    
    pl_mass = catalog["pl_masse"].to(u.kg)
    pl_rad = catalog["pl_rade"].to(u.m)
    pl_orb = catalog["pl_orbsmax"].to(u.m)
    t_sat = catalog["t_sat"].to(u.s)
    st_lum = catalog["st_lum"] * L_sun

    loss, f_xuv, _ = calculate_total_mass_loss(eta,
                                    pl_mass,
                                    R_xuv * pl_rad,
                                    pl_rad,
                                    st_lum,
                                    catalog["Lxuv/Lbol"],
                                    t_sat,
                                    (end_time * u.Gyr).to(u.s),
                                    catalog["gamma"],
                                    pl_orb)

    # Normalised mass to planet mass:
    loss = loss / pl_mass

    if output == "mass":
        colname_loss = f"Loss/pl_mass, t={end_time:.1f} Gyr"
    elif output == "fraction":
        loss = loss / protoatmosphere_mass_fraction
        colname_loss = f"Loss/{protoatmosphere_mass_fraction}protoatm., t={end_time:.1f} Gyr"
    else:
        raise ValueError("output must be 'mass' or 'fraction'")
    loss_table[colname_loss] = loss.decompose()
    loss_table[colname_loss].format = ".8f"

    colname_insol = f"insol_{end_time}_Gyr"

    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
    insolation = f_xuv / f_xuv[earth_idx]

    loss_table[colname_insol] = insolation
    loss_table[colname_insol].format = ".5f"

    return loss_table

def add_loss_colums_star_age(loss_table: Table, catalog: Table, R_xuv: float, eta: float, protoatmosphere_mass_fraction: float, output: str) -> Table:
    """Add a colum to the table with the mass loss for a planet up to the star's age
    
    Parameters
    ----------
    loss_table: Table
        the table to which the column will be added
    catalog: Table
        the original catalog with the planet data
    end_time: float
        the time at which to calculate the mass loss
    R_xuv: float
        the ratio of the XUV radius to the planet radius
    eta: float
        the heating efficiency
    protoatmosphere_mass_fraction: float
        the fraction of the planet's mass that is in the protoatmosphere
    output: str
        the type of output to calculate ("mass" or "fraction")

    Returns
    -------
    loss_table: Table
        the updated table with the added column"""
    
    pl_mass = catalog["pl_masse"].to(u.kg)
    pl_rad = catalog["pl_rade"].to(u.m)
    pl_orb = catalog["pl_orbsmax"].to(u.m)
    t_sat = catalog["t_sat"].to(u.s)
    st_lum = catalog["st_lum"] * L_sun
    st_age = catalog["st_age"].to(u.s)

    loss, f_xuv, _ = calculate_total_mass_loss(eta,
                                    pl_mass,
                                    R_xuv * pl_rad,
                                    pl_rad,
                                    st_lum,
                                    catalog["Lxuv/Lbol"],
                                    t_sat,
                                    st_age,
                                    catalog["gamma"],
                                    pl_orb)

    # Normalised mass to planet mass:
    loss = loss / pl_mass

    if output == "mass":
        colname_loss = f"Loss/pl_mass, star_age"
    elif output == "fraction":
        loss = loss / protoatmosphere_mass_fraction
        colname_loss = f"Loss/{protoatmosphere_mass_fraction}protoatm., star_age"
    else:
        raise ValueError("output must be 'mass' or 'fraction'")
    loss_table[colname_loss] = loss.decompose()
    loss_table[colname_loss].format = ".8f"

    colname_insol = f"insol_star_age"

    earth_idx = np.where(catalog["pl_name"] == "Earth")[0][0]
    insolation = f_xuv / f_xuv[earth_idx]

    loss_table[colname_insol] = insolation
    loss_table[colname_insol].format = ".5f"

    return loss_table

def main():
    catalog = nasa_exoplanet_download.main()
    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr

    output = "fraction" # "mass" or "fraction"
    R_xuv = 1.0  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01
    
    description = f"Catalog_mass_loss_for_{end_time[0]}-{end_time[-1]}_Gyr_eta-{eta}_Rxuv-{R_xuv}"


    target = "catalog" # "catalog" or "new"
    catalog = total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target, output)
    save_catalog(catalog, "AR", description, filetype="ecsv")
    
if __name__ == "__main__":
    main()