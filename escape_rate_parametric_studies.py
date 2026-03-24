from astropy.constants import L_sun,M_earth,R_earth
import numpy as np
from Plot_tools.plot_creator import Plot2D_creator, Six_3Dplot_creator, save_plot
from Calculators.function_solvers import atmospheric_escaperates_calculator
from astropy import units as u


def calculate_multiple_planets(planets, m_p, r_p, eta, l_bol, distances, r_xuv_factors, resolution):
    plot_creator = Six_3Dplot_creator(distances,r_xuv_factors, resolution)
    distances_mesh, r_xuv_factors_mesh = plot_creator.create_mesh()

    for i in range(len(planets)):
        #print(f"Calculating {planet[i]}...")
        #print(f"Mass: {m_p[i]} kg"
        #   f"\nRadius: {r_p[i]} m")
        r_xuv_mesh = r_xuv_factors_mesh * r_p[i]
        #print(f"rxuv = {r_xuv_mesh[200][200]} m")
        #print(f"Distance: {distances_mesh[200][200]} m")
        plot_creator.add_z_mesh(atmospheric_escaperates_calculator(distances_mesh,r_xuv_mesh,eta,m_p[i],r_p[i],l_bol))
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot_creator.add_norm()
    plot =plot_creator.six_window_plot(planets)
    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_r_xuv_distance_eta={eta}")

def parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,distances):
    plot_creator = Plot2D_creator(distances)
    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distances,r_xuv,eta,m_p[i],r_p[i],l_bol)
        plot_creator.append_y_axis(y_axis)
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot = plot_creator.create_2D_plot(x_label="Distance from star (AU)",y_label="Atmospheric escape rate (kg/s)",label=planets,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_standard_planets_distance_eta={eta}_rxuv_factor={r_xuv_factor}")

def parametric_r_xuv_factor(planets,m_p,r_p,r_xuv_factors,eta,l_bol,distance):
    plot_creator = Plot2D_creator(r_xuv_factors)
    for i in range(len(planets)):
        r_xuv = r_xuv_factors * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distance,r_xuv,eta,m_p[i],r_p[i],l_bol)
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"$R_{\mathrm{XUV}} / R_p$",y_label="Atmospheric escape rate (kg/s)",label=planets,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_r_xuv_factor_standard_planets_eta={eta}_distance={distance}m")

def parametric_heating_efficiency(planets,m_p,r_p,r_xuv_factor,l_bol,distance,etas):
    plot_creator = Plot2D_creator(etas)
    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distance,r_xuv,etas,m_p[i],r_p[i],l_bol)
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"Heating efficiency $\eta$",y_label="Atmospheric escape rate (kg/s)",label=planets,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_heating_efficiency_standard_planets_distance={distance}m_r_xuv_factor={r_xuv_factor}")

def parametric_l_bol(planets,m_p,r_p,r_xuv_factor,eta,l_bols,distance):
    plot_creator = Plot2D_creator(l_bols)
    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distance,r_xuv,eta,m_p[i],r_p[i],l_bols)
        plot_creator.append_y_axis(y_axis)
    plot_creator.normalize_L_sun(normalize_x_axis=True)
    plot = plot_creator.create_2D_plot(x_label=r"Bolometric Luminosity $L_{\mathrm{bol}}$ / $L_{\odot}$",y_label="Atmospheric escape rate (kg/s)",label=planets,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_l_bol_standard_planets_eta={eta}_distance={distance}m_r_xuv_factor={r_xuv_factor}")

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

    star = "Sun"
    l_bol = L_sun.value
    eta = 0.3
    distances = ((np.logspace(-1, 1, resolution))*u.AU).to_value(u.m)
    r_xuv_min = 1 ; r_xuv_max = 2   #How many times the planets atmosphere is the planets radii
    r_xuv_factors = np.logspace(np.log10(r_xuv_min), np.log10(r_xuv_max), resolution)
    calculate_multiple_planets(planets, m_p, r_p, eta, l_bol, distances, r_xuv_factors, resolution)

    r_xuv_factor = 1
    parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,distances)
    distance = (1*u.AU).to_value(u.m)
    parametric_r_xuv_factor(planets,m_p,r_p,r_xuv_factors,eta,l_bol,distance)
    etas = np.linspace(0, 1, resolution)
    parametric_heating_efficiency(planets,m_p,r_p,r_xuv_factor,l_bol,distance,etas)
    l_bols = np.logspace(-1, 1, resolution)*L_sun.value
    parametric_l_bol(planets,m_p,r_p,r_xuv_factor,eta,l_bols,distance)




if __name__=="__main__":
    main()
