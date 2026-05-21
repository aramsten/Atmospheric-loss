from astropy.constants import L_sun,M_earth,R_earth
import numpy as np
import matplotlib.pyplot as plt
from Plot_tools.plot_creator import Plot2D_creator, Six_3Dplot_creator, save_plot
from Calculators.function_solvers import atmospheric_escaperates_calculator
from astropy import units as u
from matplotlib.lines import Line2D

def calculate_multiple_planets_parametric_study_r_xuv_distance(planets, m_p, r_p, eta, l_bol, distances, r_xuv_factors, resolution):
    plot_creator = Six_3Dplot_creator(distances,r_xuv_factors, resolution,fontsize=20)
    distances_mesh, r_xuv_factors_mesh = plot_creator.create_mesh()

    for i in range(len(planets)):
        r_xuv_mesh = r_xuv_factors_mesh * r_p[i]
        plot_creator.add_z_mesh(atmospheric_escaperates_calculator(distances_mesh,r_xuv_mesh,eta,m_p[i],r_p[i],l_bol))
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot_creator.add_norm()
    plot =plot_creator.six_window_plot(planets,"Distance from star (AU)",r"$R_{\mathrm{XUV}} / R_p$",r"Atmospheric escape rate (kg/s)")
    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_r_xuv_distance_eta={eta}")

def parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,distances):
    plot_creator = Plot2D_creator(distances,fontsize=25)
    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distances,r_xuv,eta,m_p[i],r_p[i],l_bol)
        plot_creator.append_y_axis(y_axis)
    plot_creator.normalize_AU(normalize_x_axis=True)
    plot = plot_creator.create_2D_plot(x_label="Distance from star (AU)",y_label="Atmospheric escape rate (kg/s)",label=planets,x_logscale=True,y_logscale=True,view_legend=True)
    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_standard_planets_distance_eta={eta}_rxuv_factor={r_xuv_factor}")

def parametric_r_xuv_factor(planets,m_p,r_p,r_xuv_factors,eta,l_bol,distance):
    plot_creator = Plot2D_creator(r_xuv_factors,fontsize=25)
    for i in range(len(planets)):
        r_xuv = r_xuv_factors * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distance,r_xuv,eta,m_p[i],r_p[i],l_bol)
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"$R_{\mathrm{XUV}} / R_p$",y_label="Atmospheric escape rate (kg/s)",label=planets,x_logscale=True,y_logscale=True,view_legend=True)
    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_r_xuv_factor_standard_planets_eta={eta}_distance={distance}m")

def parametric_heating_efficiency(planets,m_p,r_p,r_xuv_factor,l_bol,distance,etas):
    plot_creator = Plot2D_creator(etas,fontsize=25)
    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distance,r_xuv,etas,m_p[i],r_p[i],l_bol)
        plot_creator.append_y_axis(y_axis)
    plot = plot_creator.create_2D_plot(x_label=r"Heating efficiency $\eta$",y_label="Atmospheric escape rate (kg/s)",label=planets,x_logscale=True,y_logscale=True,view_legend=True)
    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_heating_efficiency_standard_planets_distance={distance}m_r_xuv_factor={r_xuv_factor}")

def parametric_l_bol(planets,m_p,r_p,r_xuv_factor,eta,l_bols,distance):
    plot_creator = Plot2D_creator(l_bols,fontsize=25)
    for i in range(len(planets)):
        r_xuv = r_xuv_factor * r_p[i]
        y_axis = atmospheric_escaperates_calculator(distance,r_xuv,eta,m_p[i],r_p[i],l_bols)
        plot_creator.append_y_axis(y_axis)
    plot_creator.normalize_L_sun(normalize_x_axis=True)
    plot = plot_creator.create_2D_plot(x_label=r"Bolometric Luminosity $L_{\mathrm{bol}}$ / $L_{\odot}$",y_label="Atmospheric escape rate (kg/s)",label=planets,x_logscale=True,y_logscale=True,view_legend=True)
    save_plot(plot, "ST", f"atmospheric_escape_rate_parametric_study_l_bol_standard_planets_eta={eta}_distance={distance}m_r_xuv_factor={r_xuv_factor}")

def merged_1variable_studies(planets:list,m_p:list,r_p:list,r_xuv_factor:float,eta:float,l_bol:float,distance:float,normalization_factor:list):
        """Creates a plot that compined parametric studies of the distance, r_xuv_factor, heating efficiency and bolometric luminosity. The x-axis is a normalization factor that is applied to one of the varriables.
        
        Parameters
        ---
        planets: list of planet names
        m_p: list of planet masses in kg
        r_p: list of planet radii in m
        r_xuv_factor: the factor that multiplies the planet radius to get the xuv radius
        eta: heating efficiency
        l_bol: bolometric luminosity in W
        distance: distance from the star in meters
        normalization_factor: list of values to normalize the x-axis

        Output
        ---
        A plot saved to the folder "Plots"
        """
        plt.clf()
        colors = plt.cm.viridis(np.linspace(0, 1, len(planets)))
        line_dir = {
            "Distance": "-",
            r"$R_{\mathrm{XUV}} / R_p$": "--",
            "Heating efficiency": "-.",
            "Bolometric Luminosity": ":"
        }
        markers=["o","s","^","D","P","X"]
        for i in range(len(planets)):
            r_xuv = normalization_factor * r_p[i]
            y_axis = atmospheric_escaperates_calculator(distance,r_xuv,eta,m_p[i],r_p[i],l_bol)
            plt.plot(normalization_factor,y_axis,linewidth=2, color=colors[i], linestyle=line_dir[r"$R_{\mathrm{XUV}} / R_p$"], marker=markers[i], markevery=len(normalization_factor)//10+i)
            r_xuv = r_xuv_factor * r_p[i]

            y_axis = atmospheric_escaperates_calculator(distance,r_xuv,normalization_factor,m_p[i],r_p[i],l_bol)
            plt.plot(normalization_factor,y_axis,linewidth=2, color=colors[i], linestyle=line_dir["Heating efficiency"], marker=markers[i], markevery=len(normalization_factor)//10+i)

            y_axis = atmospheric_escaperates_calculator(normalization_factor,r_xuv,eta,m_p[i],r_p[i],l_bol)
            plt.plot(normalization_factor,y_axis,linewidth=2, color=colors[i], linestyle=line_dir["Distance"], marker=markers[i], markevery=len(normalization_factor)//10+i)

            y_axis = atmospheric_escaperates_calculator(distance,r_xuv,eta,m_p[i],r_p[i],normalization_factor)
            plt.plot(normalization_factor,y_axis,linewidth=2, color=colors[i], linestyle=line_dir["Bolometric Luminosity"], marker=markers[i], markevery=len(normalization_factor)//10+i)

        legend_elements = []
        for variable, linestyle in line_dir.items():
            legend_elements.append(Line2D([0], [0], color="black", lw=2, linestyle=linestyle, label=f"{variable}"))
        for i in range(len(planets)):
            color = colors[i]
            legend_elements.append(Line2D([0], [0], color=color, lw=2, label=f"{planets[i]}", marker=markers[i]))

        plt.legend(handles=legend_elements, loc="upper right", fontsize="small")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Normalization Factor",fontsize=20)
        plt.ylabel("Atmospheric escape rate (kg/s)",fontsize=20)
        plt.grid(True, which="both", ls="--", alpha=0.8)

        save_plot(plt, "ST", f"merged_1variable_study_standard_planets_eta={eta}_distance={distance}m_r_xuv_factor={r_xuv_factor}")


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
    eta = 0.1
    distances = ((np.logspace(-1, 1, resolution))*u.AU).to_value(u.m)
    r_xuv_min = 1 ; r_xuv_max = 10   #How many times the planets atmosphere is the planets radii
    r_xuv_factors = np.logspace(np.log10(r_xuv_min), np.log10(r_xuv_max), resolution)
    calculate_multiple_planets_parametric_study_r_xuv_distance(planets, m_p, r_p, eta, l_bol, distances, r_xuv_factors, resolution)

    r_xuv_factor = 1
    parametric_distance(planets,m_p,r_p,r_xuv_factor,eta,l_bol,distances)
    distance = (1*u.AU).to_value(u.m)
    parametric_r_xuv_factor(planets,m_p,r_p,r_xuv_factors,eta,l_bol,distance)
    etas = np.linspace(0, 1, resolution)
    parametric_heating_efficiency(planets,m_p,r_p,r_xuv_factor,l_bol,distance,etas)
    l_bols = np.logspace(-1, 1, resolution)*L_sun.value
    parametric_l_bol(planets,m_p,r_p,r_xuv_factor,eta,l_bols,distance)

    normalization_factor = np.logspace(np.log10(0.1), 30, resolution)
    merged_1variable_studies(planets,m_p,r_p,r_xuv_factor,eta,l_bol,distance,normalization_factor)


if __name__=="__main__":
    main()