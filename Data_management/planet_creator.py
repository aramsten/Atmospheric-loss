def make_earth_row(catalog):
    earth = catalog[:0].copy()
    earth.add_row()

    earth["pl_name"][0] = "Earth"
    earth["hostname"][0] = "Sun"
    earth["tic_id"][0] = "Sun"

    earth["pl_orbsmax"][0] = 1.0
    earth["pl_rade"][0] = 1.0
    earth["pl_masse"][0] = 1.0

    earth["st_spectype"][0] = "G2V"
    earth["st_teff"][0] = 5772
    earth["st_rad"][0] = 1.0
    earth["st_mass"][0] = 1.0
    earth["st_lum"][0] = 1.0
    earth["st_age"][0] = 4.568                  # Gyr, BouvierWadhwa2010
    earth["SpT_PM"][0] = "G2V"

    earth["Lxuv/Lbol"][0] = 0.001               # ((1.6e30*u.erg/u.s)/(L_sun.to(u.erg/u.s))).value   # ZahnleCatling2017 eq.31
    earth["Lxuv/Lbol"].format = ".5f"
    earth["t_sat"][0] = 0.1                     # Gyr, ZahnleCatling2017 under eq.31
    earth["gamma"][0] = 1.23                    # Luger2015 eq 1

    return earth
