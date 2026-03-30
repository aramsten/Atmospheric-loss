from astropy.table import vstack

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

def make_venus_row(catalog):
    venus = catalog[:0].copy()
    venus.add_row()

    venus["pl_name"][0] = "Venus"
    venus["hostname"][0] = "Sun"
    venus["tic_id"][0] = "Sun"

    venus["pl_orbsmax"][0] = 0.72 # https://www.esa.int/Science_Exploration/Space_Science/Venus_Express/Venus_compared_to_Earth
    venus["pl_rade"][0] = 0.949 # https://www.esa.int/Science_Exploration/Space_Science/Venus_Express/Venus_compared_to_Earth
    venus["pl_masse"][0] = 0.814
    """
    https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/Meet_Mercury
    https://www.esa.int/Science_Exploration/Space_Science/Venus_Express/Venus_compared_to_Earth
    """
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

    return venus

def add_planets_from_star_system(catalog, pl_names, pl_masses, pl_radii, pl_orbsmax, hostname="Sun", st_teff=5772, st_spt="G2V", st_rad=1.0, st_mass=1.0, st_lum=1.0, st_age=4.6, lq=0.001, t_sat=0.1, gamma=1.23):
    planets = catalog[:0].copy()
    for i, name in enumerate(pl_names):
        i
        planets.add_row()

        planets["pl_name"][i] = name
        planets["pl_orbsmax"][i] = pl_orbsmax[i]
        planets["pl_rade"][i] = pl_radii[i]
        planets["pl_masse"][i] = pl_masses[i]

    planets["hostname"][:] = hostname
    planets["tic_id"][:] = hostname
    planets["st_spectype"][:] = st_spt
    planets["SpT_PM"][:] = st_spt

    planets["st_teff"][:] = st_teff
    planets["st_rad"][:] = st_rad
    planets["st_mass"][:] = st_mass
    planets["st_lum"][:] = st_lum
    planets["st_age"][:] = st_age

    planets["Lxuv/Lbol"][:] = lq
    planets["Lxuv/Lbol"].format = ".5f"
    planets["t_sat"][:] = t_sat
    planets["gamma"][:] = gamma
    catalog = vstack([planets, catalog])

    return catalog