from astropy.table import vstack
from Data_management import column_creator as cc

def add_planets_from_star_system(catalog, pl_names, pl_masses, pl_radii, pl_orbsmax, hostname="Sun", st_teff=5772, st_spt="G2V", st_rad=1.0, st_mass=1.0, st_lum=1.0, st_age=4.6, lq=0.001, t_sat=0.1, gamma=1.23):
    """Add planets from a star system to the catalog with their respective properties.
    Parameters
    ---------- 
    catalog : Table
        The catalog to which the planets will be added.
    pl_names : list[str]
        A list of planet names.
    pl_masses : list[float]
        A list of planet masses in Earth masses.
    pl_radii : list[float]
        A list of planet radii in Earth radii.
    pl_orbsmax : list[float]
        A list of planet orbital semi-major axes in AU.
    hostname : str, optional
        The name of the host star (default is "Sun").
    st_teff : float, optional
        The effective temperature of the host star in Kelvin (default is 5772 K).
    st_spt : str, optional
        The spectral type of the host star (default is "G2V").
    st_rad : float, optional
        The radius of the host star in Solar radii (default is 1.0).
    st_mass : float, optional
        The mass of the host star in Solar masses (default is 1.0).
    st_lum : float, optional
        The luminosity of the host star in Solar luminosities (default is 1.0
    st_age : float, optional
        The age of the host star in Gyr (default is 4.6 Gyr
    lq : float, optional
        The ratio of X-ray to bolometric luminosity (default is 0.001).
    t_sat : float, optional
        The saturation time for the planets (default is 0.1).
    gamma : float, optional
        The gamma parameter for the planets (default is 1.23).

    Returns
    -------
    catalog : Table
        The updated catalog with the added planets.
    """
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

def add_planets_from_our_star_system(catalog,
                                     planets = ["Earth","Mercury","Venus","Mars","Moon","Jupiter","Saturn","Uranus","Neptune"]
):
    """Add planets from our star system to the catalog with their respective properties.
    
    Parameters
    ----------
    catalog : Table
        The catalog to which the planets will be added.

    Returns
    -------
    catalog : Table
        The updated catalog with the added planets.
    """

    """
    Data sources:
    https://www.esa.int/Science_Exploration/Space_Science/BepiColombo/Meet_Mercury
    https://www.esa.int/Science_Exploration/Space_Science/Venus_Express/Venus_compared_to_Earth
    https://www.esa.int/Science_Exploration/Space_Science/Mars_Express/Facts_about_Mars
    https://sci.esa.int/web/solar-system/-/35850-the-moon
    https://www.esa.int/Science_Exploration/Space_Science/Juice/Facts_about_Jupiter
    https://www.nasa.gov/wp-content/uploads/2023/05/saturn-lithograph.pdf
    https://sci.esa.int/web/solar-system/-/35653-uranus
    https://science.nasa.gov/neptune/neptune-facts/#h-size-and-distance
    """
    planet_masses = [1.0, 0.055, 0.814, 0.11, 0.0123, 318.0, 95.16, 14.536, 17.147]  # Earth masses
    planet_radii = [1.0, 0.383, 0.949, 0.533, 0.2725, 11, 9.46, 4.007, 3.887]
    planet_distance = [1.0, 0.4, 0.72, 1.52, 1.0, 5.2, 9.6, 19.182, 30]
    catalog = add_planets_from_star_system(catalog, planets, planet_masses, planet_radii, planet_distance)

    catalog = cc.add_escape_velocity(catalog, "pl_masse", "pl_rade")

    return catalog