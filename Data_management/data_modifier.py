import numpy as np

def remove_spt(catalog, colname, spt_to_remove):
    """ Removes planets with spectral type stars in string spt_to_remove, eg. "OAB" """
    spt_col = catalog[colname].astype(str)
    mask = np.ones(len(spt_col), dtype=bool)

    for spt in spt_to_remove:
        mask &= ~np.char.startswith(spt_col, spt)
    return catalog[mask]