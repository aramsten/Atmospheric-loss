from astropy.constants import G
from math import pi
import numpy as np

def luminosity_calculator(t_mesh, r_mesh):
    """Generates a mesh proposional to the luminosity of a star and normalized to our sun

    Parameters
    ---
    t_mesh: A mesh with tempatures of in Kelvin
    r_mesh: A mesh with radii normalized to our sun

    Returns 
    ---
    A mesh with luminosity normalized to our sun"""
    return r_mesh**2 * t_mesh**4

def atmospheric_escaperates_calculator(distances,r_xuv,eta,m_p,r_p,l_bol):
    """Calculates the atmospheric escape rate for a certain platet at a specified distance
    from a star. Can take inputs in the form of induvidual values or meshes.
      The output will be a mesh if the input is a mesh and a single value if the input is a single value.

    Parameters
    ---
    distances: A  distances from the star to the planet
    r_xuv: A the xuv radius
    eta: The planets heating efficiency
    m_p: The mass of the planet
    r_p: The planet radii
    l_bol: The stars Bolometric luminosity

    Returns
    ---
    Atmospheric escape rates, dimmentionless"""

    l_xuv = calculate_l_xuv(l_bol)
    f_xuv = calculate_stellar_flux(l_xuv, distances)
    m_esc_dt = eta*pi*r_xuv**2*r_p*f_xuv/(G*m_p)
    return m_esc_dt.value

def calculate_stellar_flux(l_xuv, distance):
    """Calculates the time dependant stellar flux
    
    Parameters
    ---
    l_xuv: The xuv luminosity of the star
    distance: The distance from the star

    Returns
    ---
    The stellar flux at the distance from the star"""
    f_xuv = l_xuv / (4*pi*distance**2)
    return f_xuv

def calculate_l_xuv(l_bol,t=1,t_sat=1, l_q = 1e-4):
    """Calculates the xuv-luminosity from a star
    
    Parameters
    ---
    l_bol: The bolometric luminosity of the star
    t: The age of the star
    t_sat: The saturation age of the star
    l_q: The fraction of the bolometric luminosity that is emitted in xuv

    Returns
    ---
    The xuv luminosity of the star"""
    l_q = 1e-4 #Fraction of the bolometric luminosity that is emitted in xuv.

    if t > t_sat:
        xuv_bol = l_q*(t/t_sat)**(-1.23)
    else:
        xuv_bol = l_q

    return xuv_bol*l_bol

def calculate_total_l_xuv(l_bol, l_q, t_sat, t_end, gamma):
    """Calculates the total luminosity from a star for a certain period starting at its birth.

    Parameters
    ---
    l_bol: The bolometric luminosity of the star
    l_q: The fraction of the bolometric luminosity that is emitted in xuv
    t_sat: The saturation age of the star
    t_end: The end time of the period
    gamma: The power law index for the unsaturated phase

    Returns
    ---
    The total xuv luminosity of the star over the period
    """
    l_xuv_sat = l_bol * l_q * np.minimum(t_sat, t_end)

    l_xuv_unsat = np.where(
        t_end > t_sat,
        l_bol * l_q * t_sat * (((t_end / t_sat)**(1 - gamma)) - 1) / (1 - gamma), # iff t_end > t_sat
        0.0) # iff t_end < t_sat

    return l_xuv_sat + l_xuv_unsat

def calculate_total_mass_loss(eta, m_p, r_xuv, r_p, l_bol, l_q, t_sat, t_end, gamma, distance):
    """Calculates the total mass loss from a star

    Parameters
    ---
    eta: The planets heating efficiency
    m_p: The mass of the planet
    r_xuv: A the xuv radius
    r_p: The planet radii
    l_bol: The stars Bolometric luminosity
    l_q: The fraction of the bolometric luminosity that is emitted in xuv
    t_sat: The saturation age of the star
    t_end: The end time of the period
    gamma: The power law index for the unsaturated phase
    distance: The distance from the star

    Returns
    ---
    total_mass_loss: The total mass loss from the star over the period
    total_f_xuv: The total xuv flux at the distance from the star
    total_l_xuv: The total xuv luminosity of the star over the period
    """
    total_l_xuv = calculate_total_l_xuv(l_bol, l_q, t_sat, t_end, gamma)
    total_f_xuv = calculate_stellar_flux(total_l_xuv, distance)
    total_mass_loss = eta*pi*r_xuv**2*r_p*total_f_xuv/(G*m_p)
    return total_mass_loss, total_f_xuv, total_l_xuv

def calculate_escape_velocity(m_p, r_p):
    """Calculates the escape velocity of a planet
    
    Parameters    
    ---
    m_p: The mass of the planet
    r_p: The radius of the planet
    Returns
    ---
    v_esc: The escape velocity of the planet"""
    v_esc = np.sqrt(2*G*m_p/r_p)
    return v_esc