import pyvo
import numpy as np
from Data_management import file_manager as fm
from Data_management import column_creator as cc
from Data_management import data_modifier as dm
from astropy.table import Table, vstack
from astropy import units as u
from Data_management.planet_creator import make_earth_row, add_planets_from_star_system
from astropy.constants import R_sun, L_sun, M_sun, M_earth, R_earth

def download_NASA_exoplanet_catalog(query):
  exocat = pyvo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")
  result = exocat.search(query)
  catalog = result.to_table()
  catalog["st_lum"] = 10**catalog["st_lum"]
  catalog["st_lum"].format = ".3f"
  return catalog


def apply_units(catalog):
    """ Set the correct units for the selected columns"""
    # Planet
    catalog["pl_orbsmax"] = fm.fix_unit(catalog["pl_orbsmax"], u.AU)
    catalog["pl_rade"] = fm.fix_unit(catalog["pl_rade"], u.R_earth)
    catalog["pl_masse"] = fm.fix_unit(catalog["pl_masse"], u.M_earth)

    # Star
    catalog["st_teff"] = fm.fix_unit(catalog["st_teff"], u.K)
    catalog["st_rad"] = fm.fix_unit(catalog["st_rad"], u.R_sun)
    catalog["st_mass"] = fm.fix_unit(catalog["st_mass"], u.M_sun)
    catalog["st_lum"] = fm.fix_unit(catalog["st_lum"], u.dimensionless_unscaled)
    catalog["st_age"] = fm.fix_unit(catalog["st_age"], u.Gyr)
    return catalog

def main():
    initials = "AR"

    selected_data = '''
            SELECT
              pl_name,
              hostname,
              tic_id,
              pl_orbsmax,
              pl_rade,
              pl_masse,
              st_spectype,
              st_teff,
              st_rad,
              st_mass,
              st_lum,
              st_age
            FROM
              pscomppars
            WHERE
              disc_facility LIKE '%TESS%' AND
              sy_pnum >= 1 AND
              sy_snum >= 1
            '''

    columns_to_clean = ['tic_id',
                        'pl_orbsmax',
                        'pl_rade',
                        'pl_masse',
                        'st_teff']

    spectral_types_to_remove = "OBA"

    catalog = download_NASA_exoplanet_catalog(selected_data)
    catalog = fm.filter_catalog(catalog, columns_to_clean) # Use catalog.colnames as second input if you want to clean everything, otherwise, define columns_to_clean
    catalog = apply_units(catalog)
    catalog = cc.spectral_type_from_teff(catalog, "st_teff")
    catalog = dm.remove_spt(catalog, "SpT_PM", spectral_types_to_remove)

    catalog = cc.assign_Lq_for_FGKM(catalog, "SpT_PM")
    catalog = cc.assign_saturation_for_FGKM(catalog, "SpT_PM")

    ## Used for debugging
    #for col in ["pl_rade", "pl_masse", "st_rad", "pl_orbsmax"]:
    #    print(col, catalog[col].unit, type(catalog[col].unit))
    #print(type(catalog["pl_orbsmax"]))
    #print(type(catalog["pl_orbsmax"][0]))

    catalog.sort(['tic_id'])
    #earth = make_earth_row(catalog)
    #catalog = vstack([earth, catalog])

    planets = ["Earth","Mercury","Venus","Mars"]
    """
    https://science.nasa.gov/mercury/facts/
    https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/Meet_Mercury
    
    
    https://science.nasa.gov/earth/facts/
    
    https://science.nasa.gov/mars/facts/
    https://www.esa.int/Science_Exploration/Space_Science/Mars_Express/Facts_about_Mars
    https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mars/
    
    
    """
    planet_masses = [1.0, 0.055, 0.814, 0.11]
    planet_radii = [1.0, 0.383, 0.949, 0.533]
    planet_distance = [1.0, 0.4, 0.72, 1.52]
    catalog = add_planets_from_star_system(catalog, planets, planet_masses, planet_radii, planet_distance)

    catalog = cc.add_escape_velocity(catalog, "pl_masse", "pl_rade")

    catalog = apply_units(catalog)

    return catalog

if __name__=="__main__":
    main()