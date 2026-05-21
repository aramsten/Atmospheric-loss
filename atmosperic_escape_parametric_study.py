from astropy.constants import L_sun,M_earth,R_earth
import numpy as np
from Plot_tools.plot_creator import Plot2D_creator, Six_3Dplot_creator, save_plot
from Calculators.function_solvers import calculate_total_mass_loss
from astropy import units as u

def calculate_multiple_planets_parametric_study_r_xuv_distance_eta(planets, m_p, r_p, eta, l_bol, distances, r_xuv_factors,l_q,t_sat, t_end, gamma, resolution):
    """Creates a plot of the atmospheric escape rate as a function of the distance from the star and the r_xuv_factor for multiple planets.
    
    Parameters
    ----------
    planets : list
        A list of the names of the planets.
    m_p : list
        A list of the masses of the planets in kg.
    r_p : list
        A list of the radii of the planets in m.
    eta : float
        The planets' heating efficiency.
    l_bol : float
        The stars' bolometric luminosity in W.
    distances : numpy.ndarray
        An array of distances from the star in m.
    r_xuv_factors : numpy.ndarray
        An array of factors for the XUV radius.
    l_q : float
        The fraction of the bolometric luminosity that is emitted in XUV.
    t_sat : astropy.units.Quantity
        The saturation age of the star.
    t_end : astropy.units.Quantity
        The end time of the period.
    gamma : float
        The power law index for the unsaturated phase.
    resolution : int
        The resolution of the mesh.
    """
    plot_creator = Six_3Dplot_creator(distances,r_xuv_factors, resolution,fontsize=20)
    distances_mesh, r_xuv_factors_mesh = plot_creator.create_mesh()

    for i in range(len(planets)):
        r_xuv_mesh = r_xuv_factors_mesh * r_p[i]
        mass_loss_mesh,_,_ = calculate_total_mass_loss(eta, m_p[i], r_xuv_mesh, r_p[i], l_bol, l_q, t_sat.to(u.s), t_end.to(u.s), gamma, distances_mesh)
        mass_loss_mesh = mass_loss_mesh/m_p[i]/0.01
        plot_creator.add_z_mesh(mass_loss_mesh)
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot_creator.add_norm()
    plot=plot_creator.six_window_plot(planets,"Distance from star (AU)",r"$R_{\mathrm{XUV}} / R_p$",r"$M_{loss}$ / 1%$M_\mathrm{planet}$")

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_r_xuv_distance_eta={eta}")

def parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,l_q,gamma,distances,t_end,t_sat):
    """Creates a plot of the atmospheric escape rate as a function of the distance from the star for multiple planets.
    
    parameters
    ----------
    planets : list
        A list of the names of the planets.
    m_p : list
        A list of the masses of the planets in kg.
    r_p : list
        A list of the radii of the planets in m.
    r_xuv_factor : float
        The factor for the XUV radius.
    eta : float
        The planets' heating efficiency.
    l_bol : float
        The stars' bolometric luminosity in W.
    l_q : float
        The fraction of the bolometric luminosity that is emitted in XUV.
    gamma : float
        The power law index for the unsaturated phase.
    distances : numpy.ndarray
        An array of distances from the star in m.
    t_end : astropy.units.Quantity
        The end time of the period.
    t_sat : astropy.units.Quantity
        The saturation age of the star."""
    plot_creator = Plot2D_creator(distances, fontsize=25)

    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis, _, _ = calculate_total_mass_loss(eta, m_p[i], r_xuv, r_p[i], l_bol, l_q, t_sat.to(u.s), t_end.to(u.s), gamma, distances)
        y_axis = y_axis/m_p[i]/0.01
        plot_creator.append_y_axis(y_axis)
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot = plot_creator.create_2D_plot(x_label="Distance from star (AU)",y_label=r"$M_{loss}$ / 1%$M_\mathrm{planet}$",label=planets,x_logscale=True,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_standard_planets_distance_eta={eta}_rxuv_factor={r_xuv_factor}")

def parametric_r_xuv_factor(planets,m_p,r_p,r_xuv_factors,eta,l_bol,l_q,gamma,distance,t_end,t_sat):
    """Creates a plot of the atmospheric escape rate as a function of the r_xuv_factor for multiple planets.
    
    parameters
    ----------
    planets : list
        A list of the names of the planets.
    m_p : list
        A list of the masses of the planets in kg.
    r_p : list
        A list of the radii of the planets in m.
    r_xuv_factors : numpy.ndarray
        An array of factors for the XUV radius.
    eta : float
        The planets' heating efficiency.
    l_bol : float
        The stars' bolometric luminosity in W.
    l_q : float
        The fraction of the bolometric luminosity that is emitted in XUV.
    gamma : float
        The power law index for the unsaturated phase.
    distance : numpy.ndarray
        An array of distances from the star in m.
    t_end : astropy.units.Quantity
        The end time of the period.
    t_sat : astropy.units.Quantity
        The saturation age of the star."""
    plot_creator = Plot2D_creator(r_xuv_factors, fontsize=25)

    for i in range(len(planets)):
        r_xuv = r_xuv_factors * r_p[i]
        y_axis, _, _ = calculate_total_mass_loss(eta, m_p[i], r_xuv, r_p[i], l_bol, l_q, t_sat.to(u.s), t_end.to(u.s), gamma, distance)
        y_axis = y_axis/m_p[i]/0.01
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"$R_{\mathrm{XUV}} / R_p$",y_label=r"$M_{loss}$ / 1%$M_\mathrm{planet}$",label=planets,x_logscale=True,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_standard_planets_r_xuv_factor_eta={eta}_distance={distance}m")

def parametric_heating_efficiency(planets,m_p,r_p,r_xuv_factor,etas,l_bol,l_q,gamma,distance,t_end,t_sat):
    """Creates a plot of the atmospheric escape rate as a function of the heating efficiency for multiple planets.
    
    parameters
    ----------
    planets : list
        A list of the names of the planets.
    m_p : list
        A list of the masses of the planets in kg.
    r_p : list
        A list of the radii of the planets in m.
    r_xuv_factor : float
        The factor for the XUV radius.
    etas : numpy.ndarray
        An array of heating efficiencies.
    l_bol : float
        The stars' bolometric luminosity in W.
    l_q : float
        The fraction of the bolometric luminosity that is emitted in XUV.
    gamma : float
        The power law index for the unsaturated phase.
    distance : numpy.ndarray
        An array of distances from the star in m.
    t_end : astropy.units.Quantity
        The end time of the period.
    t_sat : astropy.units.Quantity
        The saturation age of the star."""
    plot_creator = Plot2D_creator(etas, fontsize=25)

    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis, _, _ = calculate_total_mass_loss(etas, m_p[i], r_xuv, r_p[i], l_bol, l_q, t_sat.to(u.s), t_end.to(u.s), gamma, distance)
        y_axis = y_axis/m_p[i]/0.01
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"Heating efficiency $\eta$",y_label=r"$M_{loss}$ / 1%$M_\mathrm{planet}$",label=planets,x_logscale=True,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_standard_planets_heating_efficiency_r_xuv_factor={r_xuv_factor}_distance={distance}m")

def parametric_l_bol(planets,m_p,r_p,r_xuv_factor,eta,l_bols,l_q,gamma,distance,t_end,t_sat):
    """Creates a plot of the atmospheric escape rate as a function of the bolometric luminosity for multiple planets.

    parameters
    ----------
    planets : list
        A list of the names of the planets.
    m_p : list
        A list of the masses of the planets in kg.
    r_p : list
        A list of the radii of the planets in m.
    r_xuv_factor : float
        The factor for the XUV radius.
    eta : float
        The planets' heating efficiency.
    l_bols : numpy.ndarray
        An array of bolometric luminosities in W.
    l_q : float
        The fraction of the bolometric luminosity that is emitted in XUV.
    gamma : float
        The power law index for the unsaturated phase.
    distance : numpy.ndarray
        An array of distances from the star in m.
    t_end : astropy.units.Quantity
        The end time of the period.
    t_sat : astropy.units.Quantity
        The saturation age of the star."""
    plot_creator = Plot2D_creator(l_bols, fontsize=25)

    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis, _, _ = calculate_total_mass_loss(eta, m_p[i], r_xuv, r_p[i], l_bols, l_q, t_sat.to(u.s), t_end.to(u.s), gamma, distance)
        y_axis = y_axis/m_p[i]/0.01
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"Bolometric luminosity $L_{\mathrm{bol}}$",y_label=r"$M_{loss}$ / 1%$M_\mathrm{planet}$",label=planets,x_logscale=True,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_standard_planets_l_bol_r_xuv_factor={r_xuv_factor}_eta={eta}_distance={distance}m")

def parametric_starlife_our_system_study(planets,m_p,r_p,r_xuv_factor,eta,l_bol,l_q,gamma,distance,t_ends,t_sat):
    plot_creator = Plot2D_creator(t_ends/(4.603*10**9), fontsize=25)

    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis, _, _ = calculate_total_mass_loss(eta, m_p[i], r_xuv, r_p[i], l_bol, l_q, t_sat.to(u.s), t_ends.to(u.s), gamma, distance)
        y_axis = y_axis/m_p[i]/0.01
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"Time (yr)/Sun Age",y_label=r"$M_{loss}$ / 1%$M_\mathrm{planet}$",label=planets,x_logscale=True,y_logscale=True,view_legend=True)

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_standard_planets_starlife_r_xuv_factor={r_xuv_factor}_eta={eta}_distance={distance}m")

def parametric_startype_study(stars,m_p,r_p,r_xuv_factor,eta,l_bols,l_qs,gamma,distance,t_ends,t_sats):
    plot_creator = Plot2D_creator(t_ends/(4.603*10**9), fontsize=25)

    for i in range(len(stars)):
        r_xuv = r_xuv_factor * r_p
        y_axis, _, _ = calculate_total_mass_loss(eta, m_p, r_xuv, r_p, l_bols[i], l_qs[i], t_sats[i].to(u.s), t_ends.to(u.s), gamma, distance)
        y_axis = y_axis/m_p/0.01
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"Time (yr)/Sun Age",y_label=r"$M_{loss}$ / 1%$M_\mathrm{\oplus}$",label=stars,x_logscale=True,y_logscale=True,view_legend=True, markers=["o","s","^","D","P"])

    save_plot(plot, "ST", f"atmospheric_escape_parametric_study_startype_study_r_xuv_factor={r_xuv_factor}_eta={eta}_distance={distance}m_planet_Earth") 

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
    eta = 0.1 #Heating efficiency
    distances = ((np.logspace(-1, 1, resolution))*u.AU).to_value(u.m)

    star = "Sun"
    l_bol = L_sun.value
    l_q = 1e-3 #Fraction of the bolometric luminosity that is emitted in xuv.
    gamma = 1.23 #Power law index for the decay of the xuv luminos
    t_end = 4.603*u.Gyr #Time in years
    t_sat = 100*u.Myr #Saturation time in years

    stars = ["Sun (G2V)","AB Pic (K2V)","AF Lep (F8)","AU Mic (M1V)","GJ 1061 (M5.5V)"]
    t_sats = [t_sat, 100*u.Myr, 100*u.Myr, 600*u.Myr, 600*u.Myr]
    l_qs = [l_q, 1e-3, 1e-3, 1e-3, 1e-4]
    l_bols = [l_bol, 0.59387*L_sun.value, 1.94984*L_sun.value, 0.1053*L_sun.value, 0.0017*L_sun.value]

    r_xuv_min = 1 ; r_xuv_max = 10   #How many times the planets atmosphere is the planets radii
    r_xuv_factors = np.logspace(np.log10(r_xuv_min), np.log10(r_xuv_max), resolution)
    calculate_multiple_planets_parametric_study_r_xuv_distance_eta(planets, m_p, r_p, eta, l_bol, distances, r_xuv_factors,l_q,t_sat, t_end, gamma, resolution)

    r_xuv_factor = 1
    parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,l_q,gamma,distances,t_end,t_sat)
    distance = (1*u.AU).to_value(u.m)
    parametric_r_xuv_factor(planets,m_p,r_p,r_xuv_factors,eta,l_bol,l_q,gamma,distance,t_end,t_sat)
    etas = np.linspace(0, 1, resolution)
    parametric_heating_efficiency(planets,m_p,r_p,r_xuv_factor,etas,l_bol,l_q,gamma,distance,t_end,t_sat)
    l_bols = np.logspace(-1, 1, resolution)*L_sun.value
    parametric_l_bol(planets,m_p,r_p,r_xuv_factor,eta,l_bols,l_q,gamma,distance,t_end,t_sat)
    t_ends = np.logspace(np.log10(10**6), np.log10(10**11), resolution)*u.yr
    parametric_starlife_our_system_study(planets,m_p,r_p,r_xuv_factor,eta,l_bol,l_q,gamma,distance,t_ends,t_sat)

    parametric_startype_study(stars,m_p[0],r_p[0],r_xuv_factor,eta,l_bols,l_qs,gamma,distance,t_ends,t_sats)



if __name__=="__main__":
    main()
