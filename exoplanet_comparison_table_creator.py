from Data_management import nasa_exoplanet_download
from Data_management.file_manager import apply_units, save_catalog
import numpy as np
from astropy import units as u
from Data_management.planet_creator import add_planets_from_our_star_system
import exoplanet_table_creator as etc
from astropy.table import vstack

def extract_base_name(name):
    return name.split(",")[0]

def duplicate_catalog_with_variations(catalog, star, lq_values, t_sat_values, gamma_values, variation_ids, model_name):
    mask_star = catalog["hostname"] == star
    mask_sun = catalog["hostname"] == "Sun"

    star_planets = catalog[mask_star]
    sun_planets = catalog[mask_sun]

    base_catalog = vstack([sun_planets, star_planets])
    base_catalog["variation_id"] = "baseline"
    base_catalog["model_id"] = "baseline_model"

    catalog_extended = [base_catalog]


    for n in range(len(variation_ids)):
        new_tab = star_planets.copy()

        new_tab["variation_id"] = variation_ids[n]
        new_tab["model_id"] = model_name

        new_tab["Lxuv/Lbol"] = lq_values[n]
        new_tab["t_sat"] = t_sat_values[n] * u.Gyr
        new_tab["gamma"] = gamma_values[n]

        catalog_extended.append(new_tab)

    catalog = vstack(catalog_extended)
    return catalog

def change_lq_and_tsat_for_specific_star(catalog, star, lq, t_sat):
    mask = catalog["hostname"] == star
    catalog["t_sat"] = t_sat*u.Gyr
    catalog["Lxuv/Lbol"][mask] = lq
    return catalog


def main():
    initials = "AR"
    planets = ["Earth", "Mercury", "Venus", "Mars"]
    catalog = nasa_exoplanet_download.main()
    catalog = add_planets_from_our_star_system(catalog, planets)
    catalog = apply_units(catalog)

    star = "TRAPPIST-1"
    star_shortname = "T1"
    model_name = "Birky 2022"

    # Values that produce the lowest possible estimate goes first, nominal in the middle, and highest last.
    # The ordering is important as variation id's are defined separately
    suffixes = ["Birky low", "Birky", "Birky high"]
    specific_lq = [10**-(3.03+0.25), 10**-3.03, 10**-(3.03-0.23)]
    specific_t_sat = [3.14-1.46, 3.14, 3.14+2.22] # Gyr
    specific_gamma = [1.17+0.28, 1.17, 1.17-0.27]

    end_time = np.array([0.1, 0.6, 1, 5, 10])  # Gyr
    output = "fraction"  # "mass" or "fraction"
    R_xuv = 1.0  # dimensionless ratio >= 1
    eta = 0.1  # dimensionless heating efficiency
    protoatmosphere_mass_fraction = 0.01

    description = f"Catalog_comparison_{star}_for_{end_time[0]}-{end_time[-1]}_Gyr_eta-{eta}_Rxuv-{R_xuv}"

    target = "catalog"  # "catalog" or "new"

    variation_ids = ["low", "nom", "high"]

    catalog = duplicate_catalog_with_variations(catalog, star, specific_lq, specific_t_sat, specific_gamma, variation_ids, model_name)
    catalog = etc.total_loss_calculate_for_catalog(catalog, end_time, R_xuv, eta, protoatmosphere_mass_fraction, target,
                                               output)
    catalog.meta["Model"] = model_name
    catalog.meta["Compared star"] = star
    catalog.meta["Compared star shortname"] = star_shortname
    catalog.meta["suffixes"] = [suffixes]
    catalog.meta["specific_lq"] = [10 ** -(3.03 + 0.25), 10 ** -3.03, 10 ** -(3.03 - 0.23)]
    catalog.meta["specific_t_sat"] = [3.14 - 1.46, 3.14, 3.14 + 2.22]  # Gyr
    catalog.meta["specific_gamma"] = [1.17 + 0.27, 1.17, 1.17 - 0.28]

    save_catalog(catalog, initials, description, filetype="ecsv")

if __name__ == "__main__":
    main()