from astropy.table import Table
from astropy import units as u
import numpy as np
from Calculators import function_solvers as fs


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

def assign_saturation_for_FGKM(catalog, colname_spt, sat_F=0.1, sat_G=0.1, sat_K=0.1, sat_M=0.6):
    spt_col = catalog[colname_spt].astype(str)
    t_sat = np.zeros(len(spt_col))

    for i, val in enumerate(spt_col):
        if val.startswith("F"):
            t_sat[i] = sat_F # Gyr
        elif val.startswith("G"):
            t_sat[i] = sat_G # Gyr
        elif val.startswith("K"):
            t_sat[i] = sat_K # Gyr
        elif val.startswith("M"):
            t_sat[i] = sat_M # Gyr

    catalog["t_sat"] = t_sat*u.Gyr
    return catalog

def assign_gamma_for_FGKM(catalog, colname_spt, gamma_F=1.23, gamma_G=1.23, gamma_K=1.23, gamma_M=1.23):
    spt_col = catalog[colname_spt].astype(str)
    gamma = np.zeros(len(spt_col))

    for i, val in enumerate(spt_col):
        if val.startswith("F"):
            gamma[i] = gamma_F
        elif val.startswith("G"):
            gamma[i] = gamma_G
        elif val.startswith("K"):
            gamma[i] = gamma_K
        elif val.startswith("M"):
            gamma[i] = gamma_M

    catalog["gamma"] = gamma
    return catalog

def add_escape_velocity(catalog, colname_m_p, colname_r_p):
    m_p = catalog[colname_m_p].to(u.kg)
    r_p = catalog[colname_r_p].to(u.m)
    v_esc = fs.calculate_escape_velocity(m_p, r_p)
    catalog["v_esc"] = v_esc.to(u.km/u.s)
    catalog["v_esc"].format = ".1f"
    return catalog


def remove_spt(catalog, colname, spt_to_remove):
    spt_col = catalog[colname].astype(str)
    mask = np.ones(len(spt_col), dtype=bool)

    for spt in spt_to_remove:
        mask &= ~np.char.startswith(spt_col, spt)
    return catalog[mask]