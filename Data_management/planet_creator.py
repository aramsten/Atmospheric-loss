from astropy.table import vstack

def add_planets_from_star_system(catalog, pl_names, pl_masses, pl_radii, pl_orbsmax, hostname="Sun", st_teff=5772, st_spt="G2V", st_rad=1.0, st_mass=1.0, st_lum=1.0, st_age=4.6, lq=0.001, t_sat=0.1, gamma=1.23):
    planets = catalog[:0].copy()
    for i, name in enumerate(pl_names):
        planets.add_row()

        planets["pl_name"][i] = name
        planets["pl_orbsmax"][i] = pl_orbsmax[i]
        planets["pl_rade"][i] = pl_radii[i]
        planets["pl_masse"][i] = pl_masses[i]

    planets["hostname"][:] = hostname
    planets["st_spectype"][:] = st_spt
    planets["SpT_PM"][:] = st_spt

    planets["st_teff"][:] = st_teff
    planets["st_lum"][:] = st_lum
    planets["st_age"][:] = st_age

    planets["Lxuv/Lbol"][:] = lq
    planets["Lxuv/Lbol"].format = ".5f"
    planets["t_sat"][:] = t_sat
    planets["gamma"][:] = gamma
    catalog = vstack([planets, catalog])

    return catalog