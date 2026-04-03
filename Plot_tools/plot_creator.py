import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import datetime
from astropy.constants import  L_sun
from astropy import units as u

class Plot2D_creator():
    def __init__(self,x_axis):
        self.x_axis = x_axis
        self.y_axis = None
        self.error = None

    def normalize_planet_radii(self,r_p, normalize_x_axis = False, normalize_y_axis = False):
        """Normalizes the axis to planet radii"""
        if normalize_x_axis:
            self.x_axis = self.x_axis / r_p
        if normalize_y_axis:
            if isinstance(self.y_axis, list):
                self.y_axis = [y / r_p for y in self.y_axis]
            else:
                self.y_axis = self.y_axis / r_p

    def normalize_AU(self,normalize_x_axis = False, normalize_y_axis = False):
        """Normalizes the axis to AU"""
        if normalize_x_axis:
            self.x_axis = self.x_axis / (1*u.AU).to_value(u.m)
        if normalize_y_axis:
            if isinstance(self.y_axis, list):
                self.y_axis = [y / (1*u.AU).to_value(u.m) for y in self.y_axis]
            else:
                self.y_axis = self.y_axis / (1*u.AU).to_value(u.m)

    def normalize_L_sun(self,normalize_x_axis = False, normalize_y_axis = False):
        """Normalizes the axis to our sun's luminosity"""
        if normalize_x_axis:
            self.x_axis = self.x_axis / L_sun.value
        if normalize_y_axis:
            if isinstance(self.y_axis, list):
                self.y_axis = [y / L_sun.value for y in self.y_axis]
            else:
                self.y_axis = self.y_axis / L_sun.value

    def append_y_axis(self,y_axis):
        """Appends a y_axis to the object."""
        if self.y_axis is None:
            self.y_axis = [y_axis]
        else:
            self.y_axis.append(y_axis)

    def append_error(self,error):
        """Appends an error to the object."""
        if self.error is None:
            self.error = [error]
        else:
            self.error.append(error)

    def create_2D_plot(self, title = "",label="",x_label="",y_label="",x_logscale = False, y_logscale = False, plot_scatter = False, view_legend = False, view_errorbar = False):
        """Creates a plot with an x-axis, y-axis and a colorplot that is calculated from a chosen function"""
        plt.clf()
        if plot_scatter:
            plt.scatter(self.x_axis,
                self.y_axis,
                color='black',
                s = 8,
                alpha=1)
        for i in range(len(self.y_axis)):
            if not plot_scatter:
                plt.plot(self.x_axis,self.y_axis[i],label=label[i])
                if self.error is not None and view_errorbar:
                    plt.errorbar(self.x_axis,self.y_axis[i],yerr=self.error[i], fmt='o', markersize=0, capsize=3, alpha=0.5)

        if x_logscale:
            plt.xscale("log")
        if y_logscale:
            plt.yscale("log")
        if view_legend:
            plt.legend()

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid(True, which="both", ls="--", alpha=0.8)
        plt.tight_layout()

        return plt

class Plot3D_creator(Plot2D_creator):
    """A object for creating a mesh of a function. This object exist only to set up children, do not use"""

    def __init__(self,x_axis,y_axis,resolution = 400):
        super().__init__(x_axis)
        self.resolution = resolution
        self.y_axis = y_axis

    def create_mesh(self):
        """Creates a mesh of a function"""
        self.x_mesh, self.y_mesh = np.meshgrid(self.x_axis, self.y_axis)
        return self.x_mesh, self.y_mesh

    def normalize_planet_radii(self,r_p, normalize_x_axis = False, normalize_y_axis = False, normalize_z_axis = False):
        """Normalizes the radii to planet radii"""
        if normalize_x_axis:
            self.x_axis = self.x_axis / r_p
            self.x_mesh = self.x_mesh / r_p
        if normalize_y_axis:
            self.y_axis = self.y_axis / r_p
            self.y_mesh = self.y_mesh / r_p
        if normalize_z_axis:
            self.z_mesh = self.z_mesh / r_p

    def normalize_AU(self, normalize_x_axis=False, normalize_y_axis=False, normalize_z_axis=False):
        """Normalizes the axes (and mesh) to Astronomical Units."""
        # Normalize the 1D axis arrays using the base implementation
        super().normalize_AU(normalize_x_axis, normalize_y_axis)

        # Also normalize mesh coordinates for 3D plots (pcolormesh uses mesh grids)
        if normalize_x_axis:
            self.x_mesh = self.x_mesh / (1*u.AU).to_value(u.m)
        if normalize_y_axis:
            self.y_mesh = self.y_mesh / (1*u.AU).to_value(u.m)
        if normalize_z_axis:
            self.z_mesh = self.z_mesh / (1*u.AU).to_value(u.m)

    def set_z_mesh(self,z_mesh):
        """Sets the z_mesh to a calculated mesh of a function"""
        self.z_mesh = z_mesh

    def add_norm(self, vmin=None, vmax=None):
        """Adds a LogNorm normalized to the closest powers of ten for all added meshes."""
        if vmin is None:
            vmin = 10**np.floor(np.log10(np.min(self.z_mesh)))
        if vmax is None:
            vmax = 10**np.ceil(np.log10(np.max(self.z_mesh)))

        self.norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    def create_3D_plot(self, title = "",label="",x_label="",y_label="",color_label="", x_logscale = False, y_logscale = False, plot_scatter = False, plt_colorbar = True):
        """Creates a plot with an x-axis, y-axis and a colorplot that is calculated from a chosen function"""
        plt.clf()
        if plot_scatter:
            plt.scatter(self.x_axis,
                self.y_axis,
                color='black',
                s = 8,
                alpha=1)

        mesh = plt.pcolormesh(
            self.x_mesh, self.y_mesh, self.z_mesh,
            shading='auto',
            cmap='viridis',
            #cmap='magma',
            norm=self.norm, # Normalises colormap to max/min mesh Luminosity
            alpha=0.6
        )
        if x_logscale:
            plt.xscale("log")
        if y_logscale:
            plt.yscale("log")

        plt.xlim(np.min(self.x_axis), np.max(self.x_axis))
        plt.ylim(np.min(self.y_axis), np.max(self.y_axis))

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        if plt_colorbar == True:
            plt.colorbar(mesh, label=color_label)
        plt.grid(True, which="both", ls="--", alpha=0.8)
        plt.tight_layout()

        return plt

class Six_3Dplot_creator(Plot3D_creator):

    def __init__(self,x_axis,y_axis,resolution = 400):
        super().__init__(x_axis,y_axis,resolution)
        self.z_meshes = []

    def add_z_mesh(self,z_mesh):
        """Adds a z_mesh to the list of z_meshes"""
        self.z_meshes.append(z_mesh)

    def add_norm(self, vmin=None, vmax=None):
        """Adds a LogNorm normalized to the closest powers of ten for all added meshes."""
        if vmin is None:
            vmin = 10**np.floor(np.log10(np.min(self.z_meshes)))
        if vmax is None:
            vmax = 10**np.ceil(np.log10(np.max(self.z_meshes)))

        self.norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    def six_window_plot(self,planet):
        fig, axs = plt.subplots(2, 3, figsize=(12, 8), constrained_layout=True)

        for ax, data, title in zip(axs.flat, self.z_meshes, planet):
            mesh = ax.pcolormesh(
                self.x_mesh, self.y_mesh, data,
                shading='auto',
                cmap='viridis',
                norm=self.norm,
                alpha=0.6
            )
            ax.set_title(title)
            ax.set_xlabel("Distance to star (AU)")
            ax.set_ylabel(r"$R_{\mathrm{XUV}} / R_p$")
            ax.set_xscale("log")
            ax.grid(True, which="both", ls="--", alpha=0.8)

        fig.colorbar(mesh, ax=axs, orientation='vertical', fraction=0.02, pad=0.02, label="Atmospheric Escape Rate (kg/s)")
        return plt



def save_plot(plt,initials,plot_name, folder="Plots"):
  """Saves a plot on the computer"""
  today = datetime.datetime.today()
  today = today.strftime("%y%m%d_%H.%M")

  # Creates folder if not existing:
  os.makedirs(folder, exist_ok=True)

  plt.savefig(f"{folder}/{today}_{initials}_Plot_{plot_name}.pdf", bbox_inches='tight')