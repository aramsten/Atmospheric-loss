from astroquery.mast import Catalogs
import numpy as np
import os

def download_TIC_catalog(Teff,rad,Tmag):
  """Dowloads stars data from the TIC-Catalog. Note: If the set of data from the
  Teff: An array with the lower and upper bound of tempatures (int) of the stars that are downloaded. Ex [5000, 6300]
  rad: An array with the lower and upper bound of radii (int) of the stars that are downloaded. Ex [0.9, 1.1]
  Tmag:
  Output:
     """
  catalog = Catalogs.query_criteria(
      catalog="TIC",
      Teff=Teff,
      rad=rad,
      Tmag=Tmag,
      objType="STAR"
  )
  return catalog

def filter_catalog(catalog):
  """Filters the catalog so that only the coloms with ID, Teff, rad and lum is left
  Catalog: A table with TIC-data
  Output: A filtered table with TIC-data"""
  catalog = catalog[['ID',
                    'Teff',
                    'rad',
                    'lum']]

  mask_rad = (~catalog['rad'].mask) & (~np.isnan(catalog['rad']))
  mask_teff = (~catalog['Teff'].mask) & (~np.isnan(catalog['Teff']))
  mask_lum = (~catalog['lum'].mask) & (~np.isnan(catalog['lum']))
  mask_id = (~catalog['ID'].mask)

  mask_all = mask_rad & mask_teff & mask_lum & mask_id
  catalog = catalog[mask_all]

  return catalog

def save_catalog(catalog):
  """Saves a catalog on the computer
  Catalog: A table with TIC-data
  """
  catalog.write('catalog.csv', format='csv', overwrite = True)


def main():
    Teff=[5950, 6000]
    rad=[0.9, 1.1]
    Tmag=[6, 12]

    catalog = download_TIC_catalog(Teff,rad,Tmag)
    catalog = filter_catalog(catalog)
    catalog.sort(['ID'])
    save_catalog(catalog)


if __name__=="__main__":
    main()