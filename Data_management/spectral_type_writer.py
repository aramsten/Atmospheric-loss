from astropy.table import Table
from astropy import units as u
import numpy as np

def spectral_type_from_teff(catalog, teff_colname):
    pm_table = Table.read("Data_management/External_data/EEM_dwarf_UBVIJHK_colors_Teff.txt",
                          format='ascii',
                          header_start = 23)
    pm_table = pm_table["SpT", "Teff"]

    teff_col = catalog[teff_colname]
    teff_array = np.array(teff_col, dtype=float)
    diff = np.abs(pm_table["Teff"][:, None] - teff_array[None, :])
    idx = np.argmin(diff, axis=0)

    catalog["SpT_PM"] = pm_table["SpT"][idx]
    return catalog

def assign_Lq_for_FGKM(catalog, colname_spt):
    spt_col = catalog[colname_spt].astype(str)
    M_mask = np.char.startswith(spt_col, "M")
    M_subtype = np.full(len(spt_col), np.nan)

    for i, val in enumerate(spt_col):
        if val.startswith("M"):
            num = ""
            for ch in val[1:]:
                if ch.isdigit() or ch == ".":
                    num += ch
                else:
                    break
            M_subtype[i] = float(num)

    Lq_fgk_M00_M40 = 1e-3
    Lq_M45_M100 = 1e-4

    ratio = np.where(
        (M_mask) & (M_subtype > 4.0),
        Lq_M45_M100,
        Lq_fgk_M00_M40)

    catalog["Lxuv/Lbol"] = ratio
    return catalog

def assign_saturation_for_FGKM(catalog, colname_spt):
    spt_col = catalog[colname_spt].astype(str)

    t_sat = np.zeros(len(spt_col))
    gamma = np.zeros(len(spt_col))

    for i, val in enumerate(spt_col):
        if val.startswith("F"):
            t_sat[i] = 0.1 # Gyr
            gamma[i] = 1.23
        elif val.startswith("G"):
            t_sat[i] = 0.1 # Gyr
            gamma[i] = 1.23
        elif val.startswith("K"):
            t_sat[i] = 0.1 # Gyr
            gamma[i] = 1.23
        elif val.startswith("M"):
            t_sat[i] = 0.6 # Gyr
            gamma[i] = 1.23

    catalog["t_sat"] = t_sat*u.Gyr
    catalog["gamma"] = gamma
    return catalog

def remove_spt(catalog, colname, spt_to_remove):
    spt_col = catalog[colname].astype(str)
    mask = np.ones(len(spt_col), dtype=bool)

    for spt in spt_to_remove:
        mask &= ~np.char.startswith(spt_col, spt)
    return catalog[mask]