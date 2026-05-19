from astropy.constants import L_sun,M_earth,R_earth
import numpy as np
from Plot_tools.plot_creator import Plot3D_creator,Plot2D_creator, save_plot
from Calculators.function_solvers import calculate_total_mass_loss
from astropy import units as u

def parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,l_q,gamma,distances,t_end,t_sat):
    plot_creator = Plot2D_creator(distances)

    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = np.array([])
        for distance in distances:
            y_value, _, _ = calculate_total_mass_loss(eta, m_p[i], r_xuv, r_p[i], l_bol, l_q, t_sat.to(u.s), t_end.to(u.s), gamma, distance)
            y_axis = np.append(y_axis, y_value.value)
        y_axis = y_axis/m_p[i]/0.001
        plot_creator.append_y_axis(y_axis)
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot = plot_creator.create_2D_plot(x_label="Distance from star (AU)",y_label=r"$M_{loss}$ / 1%$M_{\oplus}$",label=planets,x_logscale=True,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_standard_planets_distance_eta={eta}_rxuv_factor={r_xuv_factor}")


def main():
    resolution = 400

    planets = ["Mercury","Venus","Sub Earth","Earth","Super Earth","Mars"]
    m_mercury = 0.055*M_earth.value
    m_venus = 4.87*10**24
    m_earth = M_earth.value
    m_super_earth = 5*M_earth.value
    m_sub_earth = 0.1*M_earth.value
    m_mars = 0.11*M_earth.value
    m_p = np.array([m_mercury, m_venus, m_sub_earth, m_earth, m_super_earth, m_mars])

    r_mercury = 0.38*R_earth.value
    r_venus = 6.052*10**6
    r_sub_earth = 0.46*R_earth.value
    r_earth = R_earth.value
    r_super_earth = 2.71*R_earth.value
    r_mars = 0.53*R_earth.value
    r_p = np.array([r_mercury, r_venus, r_sub_earth, r_earth, r_super_earth, r_mars])

    r_xuv_factor = 1
    eta = 0.3 #Heating efficiency
    distances = ((np.logspace(-1, 1, resolution))*u.AU).to_value(u.m)

    star = "Sun"
    l_bol = L_sun.value
    l_q = 1e-3 #Fraction of the bolometric luminosity that is emitted in xuv.
    gamma = 1.23 #Power law index for the decay of the xuv luminos
    t_end = 4.603*u.Gyr #Time in years
    t_sat = 100*u.Myr #Saturation time in years

    parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,l_q,gamma,distances,t_end,t_sat)

if __name__=="__main__":
    main()
