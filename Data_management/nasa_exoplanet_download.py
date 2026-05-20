import pyvo
import numpy as np
from Data_management import file_manager as fm
from Data_management import column_creator as cc
from astropy.table import Table, vstack
from Data_management.planet_creator import add_planets_from_our_star_system
from astropy.constants import R_sun, L_sun, M_sun, M_earth, R_earth

def download_NASA_exoplanet_catalog(query):
  """Downloads the NASA exoplanet catalog form the exoplanet archive using the provided query and returns it as an astropy Table.
  
  Parameters
  ----------  
  query : str
      The SQL query to be executed on the NASA exoplanet archive.
  
  Returns
  -------
  catalog : Table
      An astropy Table containing the results of the query.
  """
  exocat = pyvo.dal.TAPService("https://exoplanetarchive.ipac.caltech.edu/TAP")
  result = exocat.search(query)
  catalog = result.to_table()
  catalog["st_lum"] = 10**catalog["st_lum"]
  catalog["st_lum"].format = ".5f"
  return catalog

def main():
    initials = "AR"

    # To only use planets confirmed by TESS, use --> disc_facility LIKE '%TESS%' AND <-- under WHERE
    selected_data = '''
            SELECT
              pl_name,
              hostname,
              pl_orbsmax,
              pl_rade,
              pl_masse,
              st_spectype,
              st_teff,
              st_lum,
              st_age
            FROM
              pscomppars
            WHERE
              sy_pnum >= 1 AND
              sy_snum = 1
            '''

    columns_to_clean = ['pl_orbsmax',
                        'pl_rade',
                        'pl_masse',
                        'st_teff',
                        'st_lum',
                        'st_age'
                        ]

    spectral_types_to_remove = "OBAT"

    catalog = download_NASA_exoplanet_catalog(selected_data)
    catalog = fm.filter_catalog(catalog, columns_to_clean) # Use catalog.colnames as second input if you want to clean everything, otherwise, define columns_to_clean
    catalog = fm.apply_units(catalog)
    catalog = cc.spectral_type_from_teff(catalog, "st_teff")
    catalog = fm.remove_spt(catalog, "SpT_PM", spectral_types_to_remove)

    catalog = cc.assign_Lq_for_FGKM(catalog, "SpT_PM")
    catalog = cc.assign_saturation_for_FGKM(catalog, "SpT_PM")
    catalog = cc.assign_gamma_for_FGKM(catalog, "SpT_PM")

    catalog.sort(['pl_name'])

    return catalog

if __name__=="__main__":
    main()