import os
import datetime
import numpy as np
from astropy import units as u

def clean_column(table, colname):
    col = table[colname]

    if col.dtype.kind in ['U', 'S', 'O']: # Checks for alphabetic rows, U=Unicode, S=Byte-string, O=Object
        col_str = col.astype(str)         # Creates string from columns
        col_str = np.char.strip(col_str)  # Removes spaces from start and end of string
        mask_not_empty = (col_str != '')  # Creates mask that removes empty rows, necessary for rows containing only ''
    else:
        mask_not_empty = ~np.isnan(col)   # Creates mask that removes rows that is Nan, "not a number"
    mask_valid = ~col.mask & mask_not_empty
    return mask_valid

def filter_catalog(cat, columns_to_clean):
    mask = np.ones(len(cat), dtype=bool)  # Creates a boolean mask that is initially true for all rows
    for colname in columns_to_clean:
        mask &= clean_column(cat, colname) # Updates mask for every colname, only keeps columns containg info
    clean_cat = cat[mask]
    return clean_cat

def remove_spt(catalog, colname, spt_to_remove):
    """ Removes planets with spectral type stars in string spt_to_remove, eg. "OAB" """
    spt_col = catalog[colname].astype(str)
    mask = np.ones(len(spt_col), dtype=bool)

    for spt in spt_to_remove:
        mask &= ~np.char.startswith(spt_col, spt)
    return catalog[mask]

def fix_unit(col, unit):
    # If it is a Quantity
    if hasattr(col, "unit") and col.unit is not None:
        try:
            #print(f"1: {col.name}")
            return col.to(unit) # Used for debug
        except Exception:
            # Broken Quantity (UnrecognizedUnit etc.)
            #print(f"2: {col.name}") # Used for debug
            return col.value * unit
    # No unit assigned
    return col * unit

def apply_units(catalog):
    """ Set the correct units for the selected columns"""
    # Planet
    catalog["pl_orbsmax"] = fix_unit(catalog["pl_orbsmax"], u.AU)
    catalog["pl_rade"] = fix_unit(catalog["pl_rade"], u.R_earth)
    catalog["pl_masse"] = fix_unit(catalog["pl_masse"], u.M_earth)

    # Star
    catalog["st_teff"] = fix_unit(catalog["st_teff"], u.K)
    catalog["st_lum"] = fix_unit(catalog["st_lum"], u.dimensionless_unscaled)
    catalog["st_age"] = fix_unit(catalog["st_age"], u.Gyr)
    return catalog

def save_catalog(catalog, initials, description="Catalog_NASA_exoplanet",folder="Tables", filetype="csv"):
    today = datetime.datetime.today()
    today = today.strftime("%y%m%d_%H.%M")

    # Creates folder if not existing:
    os.makedirs(folder, exist_ok=True)

    filename = f"{today}_{initials}_{description}.{filetype}"
    catalog.write(f"{folder}/{filename}", format=filetype, overwrite=True)