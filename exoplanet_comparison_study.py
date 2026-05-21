from astropy.constants import L_sun, M_earth, R_earth, G
import matplotlib.pyplot as plt
from Plot_tools.plot_creator import save_plot
import numpy as np
from astropy.io import ascii
from astropy import units as u
from Calculators import function_solvers as fs
from escape_velocity_mass_loss_studies import draw_shoreline

def atmosphere_color(name):
    """
    Returns color based on atmospheric status of known planets
    green <==> thick
    orange <==> thin
    red <==> airless
    gray <==> unknown
    """
    name = name.replace('"', '')

    # Solar system
    if name == "Earth": return "green"
    if name == "Venus": return "green"
    if name == "Mars": return "orange"
    if name == "Jupiter": return "green"
    if name == "Moon": return "red"
    if name == "Mercury": return "red"
    if name == "Saturn": return "green"
    if name == "Uranus": return "green"
    if name == "Neptune": return "green"

    # TRAPPIST-1
    if name == "TRAPPIST-1 b": return "red"
    if name == "TRAPPIST-1 c": return "red"
    if name == "TRAPPIST-1 d": return "red"

    return "gray"

def plot_cosmic_shoreline_regular(catalog, colname, x_label, y_label, eta, model_name, star, star_shortname, fontsize):
    base = catalog[catalog["variation_id"] == "baseline"]
    model = catalog[catalog["model_id"] == model_name]

    x_base = base["v_esc"]
    y_base = base[colname]

    fig, ax = plt.subplots()

    # Baseline
    for i, name in enumerate(base["pl_name"]):
        color = atmosphere_color(name)
        ax.scatter(x_base.value[i], y_base[i],
                   color=color, marker="o", s=40, zorder=3)

        ax.text(x_base.value[i] * 1.02,
                y_base[i] * 0.8,
                short_label(name, star, star_shortname),
                fontsize=7,
                color="black",
                fontweight="bold",
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0))

    # Compared model
    for name in base["pl_name"]:
        rows = model[model["pl_name"] == name]
        if len(rows) == 0:
            continue

        y_vals = np.array([r[colname] for r in rows])
        y_min, y_max = np.min(y_vals), np.max(y_vals)

        nom_mask = rows["variation_id"] == "nom"
        nom = rows[nom_mask][0]
        y_nom = nom[colname]

        x_nom = nom["v_esc"]

        color = atmosphere_color(name)

        ax.errorbar(
            x_nom, y_nom,
            yerr=[[y_nom - y_min], [y_max - y_nom]],
            fmt='s',
            color=color,
            markersize=6,
            capsize=3,
            zorder=4
        )

        ax.text(x_nom * 1.02,
                y_nom * 1.05,
                short_label(name, star, star_shortname),
                fontsize=7,
                color="black",
                fontweight="bold",
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5))

    # Axis settings
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label,fontsize=fontsize)
    ax.set_ylabel(y_label,fontsize=fontsize)

    ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.7)

    plot_margin = 0.07
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    ax.set_xlim(x_min * (1-plot_margin) , x_max * (1+plot_margin))
    ax.set_ylim(y_min * (1-plot_margin), y_max * (1+plot_margin))

    # SHoreline
    shoreline_position_text = draw_shoreline(ax, catalog, x_base, y_base, color="blue", alpha=0.5)
    x_pos = 4.2
    y_pos = 0.23
    ax.text(x_pos, y_pos,
            r"$I_\mathrm{XUV} \propto v_\mathrm{esc}^4$",
            fontsize=10,
            va="center",
            ha="center",
            rotation=24,
            zorder=2,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.9, pad=1.2)
            )

    # Legends
    add_legends(ax, model_name)

    return fig, shoreline_position_text


def plot_cosmic_shoreline_models(catalog, colname, x_label, y_label, eta, model_name, star, star_shortname,fontsize):
    base = catalog[catalog["variation_id"] == "baseline"]
    model = catalog[catalog["model_id"] == model_name]

    earth_idx = np.where(base["pl_name"] == "Earth")[0][0]
    v_earth = base["v_esc"].quantity[earth_idx]

    mass_base = base["pl_masse"].to(u.kg)
    rad_base  = base["pl_rade"].to(u.m)
    rho_base  = 3 * mass_base / (4 * np.pi * rad_base**3)
    rho_earth = rho_base[earth_idx]

    x_base = ((base["v_esc"] / v_earth)**3 * (rho_base / rho_earth)**0.5).decompose()
    y_base = base[colname]

    fig, ax = plt.subplots()

    # Baseline
    for i, name in enumerate(base["pl_name"]):
        color = atmosphere_color(name)
        ax.scatter(x_base.value[i], y_base[i],
                   color=color, marker="o", s=40, zorder=3)

        ax.text(x_base.value[i] * 1.05,
                y_base[i] * 0.8,
                short_label(name, star, star_shortname),
                fontsize=7,
                color="black",
                fontweight="bold",
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0))

    # Compared model
    for name in base["pl_name"]:
        rows = model[model["pl_name"] == name]
        if len(rows) == 0:
            continue

        y_vals = np.array([r[colname] for r in rows])
        y_min, y_max = np.min(y_vals), np.max(y_vals)

        nom_mask = rows["variation_id"] == "nom"
        nom = rows[nom_mask][0]
        y_nom = nom[colname]

        v_nom = nom["v_esc"]* catalog["v_esc"].unit
        rho_nom = 3 * (nom["pl_masse"] * M_earth) / (4 * np.pi * (nom["pl_rade"] * R_earth)**3)
        x_nom = ((v_nom / v_earth)**3 * (rho_nom / rho_earth)**0.5).decompose()

        color = atmosphere_color(name)

        ax.errorbar(
            x_nom.value, y_nom,
            yerr=[[y_nom - y_min], [y_max - y_nom]],
            fmt='s',
            color=color,
            markersize=6,
            capsize=3,
            zorder=4
        )

        ax.text(x_nom.value * 1.05,
                y_nom * 1.05,
                short_label(name, star, star_shortname),
                fontsize=7,
                color="black",
                fontweight="bold",
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5))

    # Axis settings
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label,fontsize=fontsize)
    ax.set_ylabel(y_label,fontsize=fontsize)

    ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.7)

    plot_margin = 0.07
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    ax.set_xlim(x_min * (1-plot_margin) , x_max * (1+plot_margin))
    ax.set_ylim(y_min * (1-plot_margin), y_max * (1+plot_margin))

    # Mass loss lines
    mars_loss = (catalog[catalog['pl_name'] == 'Mars']['Loss/0.01protoatm., star_age'][0])*1e-2
    mass_fractions = [2e-4, 2e-3, 1e-2, 2e-1, mars_loss]
    draw_mass_loss_lines_for_comparison(ax, base, eta, x_base, y_base, mass_fractions, earth_idx)

    # Legends
    add_legends(ax, model_name)

    return fig

def draw_mass_loss_lines_for_comparison(ax, catalog, eta, x_axis, y_axis, fractions, earth_idx):
    mars_loss = (catalog[catalog['pl_name'] == 'Mars']['Loss/0.01protoatm., star_age'][0])*1e-2

    d = catalog["pl_orbsmax"].to(u.m)
    lum = catalog["st_lum"]*L_sun
    lq = catalog["Lxuv/Lbol"]
    t_sat = catalog["t_sat"].to(u.s)
    age = catalog["st_age"].to(u.s)
    gamma = catalog["gamma"]

    d_earth = d[earth_idx]
    lum_sun = lum[earth_idx]
    lq_sun = lq[earth_idx]
    t_sat_sun = t_sat[earth_idx]
    age_sun = age[earth_idx]
    gamma_sun = gamma[earth_idx]

    E_xuv_sun = fs.calculate_total_l_xuv(lum_sun, lq_sun, t_sat_sun, age_sun, gamma_sun)

    v_esc = catalog["v_esc"].to(u.m/u.s)
    rho = 3*catalog["pl_masse"].to(u.kg)/(4*np.pi*(catalog["pl_rade"].to(u.m))**3)
    v_earth = v_esc[earth_idx]
    rho_earth = rho[earth_idx]

    xmin, xmax = ax.get_xlim()
    x_vals = np.logspace(np.log10(xmin), np.log10(xmax), 200)

    ax = plt.gca()

    for f in fractions:
        m = (d_earth ** 2 / (eta * E_xuv_sun)) * (8 * np.pi / (3 * G))**0.5 * f *(v_earth) ** 3 * (rho_earth)**0.5
        y_line = m * x_vals
        ax.plot(x_vals, y_line, linestyle="--", color="blue", alpha=0.5, zorder=1)
        x_pos = x_vals[15] if f != 0.0002 else x_vals[80]
        y_pos = y_line[15] if f != 0.0002 else y_line[82]
        ax.text(x_pos, y_pos,
                r"$\Delta M_\mathrm{Mars}$" if f == mars_loss else rf"$\Delta M = {f}$",
                fontsize=8,
                va="center",
                ha="center",
                rotation =18,
                zorder=2,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1)
                )

def add_legends(ax, model_name):
    # Legend 1: Marker shape (baseline vs model)
    baseline_proxy = plt.Line2D([], [], marker="o", color="black", linestyle="None", label="Baseline")
    model_proxy = plt.Line2D([], [], marker="s", color="black", linestyle="None", label=f"{model_name}")
    legend1 = ax.legend(handles=[baseline_proxy, model_proxy],loc="lower left", fontsize=8, frameon=True)
    ax.add_artist(legend1)

    # Legend 2: Colour (atmospheric status)
    atm_green = plt.Line2D([], [], marker="o", color="green", linestyle="None", label="Thick atmosphere")
    atm_orange = plt.Line2D([], [], marker="o", color="orange",linestyle="None", label="Thin atmosphere")
    atm_red = plt.Line2D([], [], marker="o", color="red", linestyle="None", label="Airless")
    atm_gray = plt.Line2D([], [], marker="o", color="gray", linestyle="None", label="Unknown")
    ax.legend(handles=[atm_green, atm_orange, atm_red, atm_gray], loc="upper left", fontsize=8, frameon=True)


def short_label(name, star, star_shortname):
    """
    Returnerar kortare label för plot, ändrar inte catalog
    exempelvis blir TRAPPIST-1 e omvandlat till T1-e
    """
    if name.startswith(f"{star} "):
        suffix = name.split(f"{star} ")[1]
        return f"{star_shortname}-{suffix}"
    else:
        return name

def main():
    initials = "AR"

    table_name = "260521_23.17_AR_Catalog_comparison_TRAPPIST-1_for_0.1-10.0_Gyr_eta-0.1_Rxuv-1.0.ecsv"
    catalog = ascii.read(f"Tables/{table_name}")
    fontsize = 16
    eta = catalog.meta["eta"]
    R_xuv = catalog.meta["R_xuv"]
    model_name = catalog.meta["Model"]
    star = catalog.meta["Compared star"]
    star_shortname = catalog.meta["Compared star shortname"]

    plot_name = f"cosmic_shoreline_comparison_TRAPPIST-1_at_stars_age_lost-rxuv={R_xuv}-eta={eta}"

    plot = plot_cosmic_shoreline_models(
        catalog,
        "insol_star_age",
        r"$x = (v_\mathrm{esc}/v_\mathrm{esc\oplus})^3 (\rho / \rho_\oplus)^{0.5}$",
        "Insolation relative to Earth",
        eta,
        model_name,
        star,
        star_shortname,
        fontsize
    )
    save_plot(plot, initials, plot_name)

    plot, shoreline_position_text = plot_cosmic_shoreline_regular(
        catalog,
        "insol_star_age",
        "Escape velocity (km/s)",
        "Insolation relative to Earth",
        eta,
        model_name,
        star,
        star_shortname,
        fontsize
    )
    save_plot(plot, initials, f"{shoreline_position_text}_{plot_name}")

if __name__ == "__main__":
    main()