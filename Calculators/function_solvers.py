from astropy.constants import G
from astropy import units as u ## Used for some debugging functions, may be removed if debug is deleted
from math import pi
import numpy as np

def luminosity_calculator(t_mesh, r_mesh):
    """Generates a mesh proposional to the luminosity of a star and normalized to our sun
    t_mesh: A mesh with tempatures of in Kelvin
    r_mesh: A mesh with radii normalized to our sun
    Output: A mesh with luminosity normalized to our sun"""
    return r_mesh**2 * t_mesh**4

def atmospheric_escaperates_calculator(distances,r_xuv,eta,m_p,r_p,l_bol):
    """Calculates the atmospheric escape rate for a certain platet at a specified distance
    from a star.
     r_xuv_factors: A mesh with the dimentionless factors of the planets radii that the xuv-radiation reaches. Ex 1.5 means that the xuv-radiation reaches 1.5 times the planets radii."""
    l_xuv = calculate_l_xuv(l_bol)
    f_xuv = calculate_stellar_flux(l_xuv, distances)
    m_esc_dt = eta*pi*r_xuv**2*r_p*f_xuv/(G*m_p)
    return m_esc_dt.value

def calculate_stellar_flux(l_xuv, distance):
    """Calculates the time dependant stellar flux"""
    f_xuv = l_xuv / (4*pi*distance**2)
    return f_xuv

def calculate_l_xuv(l_bol,t=1,t_sat=1):
    """Calculates the xuv-luminosity from a star"""
    l_q = 1e-4 #Fraction of the bolometric luminosity that is emitted in xuv.

    if t > t_sat:
        xuv_bol = l_q*(t/t_sat)**(-1.23)
    else:
        xuv_bol = l_q

    return xuv_bol*l_bol

def calculate_total_l_xuv(l_bol, l_q, t_sat, t_end, gamma):
    """Calculates the total luminosity from a star. Use SI-units."""
    l_xuv_sat = l_bol * l_q * np.minimum(t_sat, t_end)

    l_xuv_unsat = np.where(
        t_end > t_sat,
        l_bol * l_q * t_sat * (((t_end / t_sat)**(1 - gamma)) - 1) / (1 - gamma), # iff t_end > t_sat
        0.0) # iff t_end < t_sat

    return l_xuv_sat + l_xuv_unsat

def calculate_total_mass_loss(eta, m_p, r_xuv, r_p, l_bol, l_q, t_sat, t_end, gamma, distance):
    """Calculates the total mass loss from a star"""
    total_l_xuv = calculate_total_l_xuv(l_bol, l_q, t_sat, t_end, gamma)
    total_f_xuv = calculate_stellar_flux(total_l_xuv, distance)
    total_mass_loss = eta*pi*r_xuv**2*r_p*total_f_xuv/(G*m_p)
    #print(f"unit for mass loss at {t_end.to(u.Gyr)}: {total_mass_loss.decompose().unit}") ## Used for debugging
    return total_mass_loss, total_f_xuv, total_l_xuv

def calculate_escape_velocity(m_p, r_p):
    """Calculates the escape velocity"""
    v_esc = np.sqrt(2*G*m_p/r_p)
    return v_esc