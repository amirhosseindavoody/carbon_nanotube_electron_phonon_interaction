# this program reads and plots the calculated exciton scattering rates
# in CNT conduction band.
import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools

# prepare for ploting
figure_list = list() # this will hold the list of all the figures
plt.ion() # enables interactive mode for ploting. No plt.show() is needed!
colors = itertools.cycle(["r", "b", "g", "k", "c", "m", "y"]) # iterative list for cycling through colors when ploting. use color=next(colors) to get the next color on the list.
styles = itertools.cycle(["solid", "dashed", "dash_dot"]) # iterative list for cycling through line styles when ploting. use linestyle=next(styles) to get the next color on the list.

# set some constants
eV = 1.6e-19 #[Jouls]
eV_to_inverse_cm = 8065.54429
hbar = 1.054e-34 #[Jouls.second]
vF = 1.0e6

directory = "/home/amirhossein/research/exciton/data/transfer_rates/final_result/"
cnt_name = "cnt1"

################################################################################
# load the calculated scattering rates for emission and absorption process
data = np.loadtxt(directory+cnt_name+".exciton_phonon_scattering_rate_emission.dat", skiprows=0)
energy_mesh_emission = data[0,:]/eV
scattering_rate_emission = data[1,:].T

energy_mesh_emission = energy_mesh_emission - np.amin(energy_mesh_emission)

################################################################################
# plot the electron-phonon scattering rates

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(1,1,1)

axes.plot(energy_mesh_emission, scattering_rate_emission[:], color="black", linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)

axes.set_yscale('log')

xmin = min(energy_mesh_emission)
xmax = max(energy_mesh_emission)
axes.set_xlim([xmin,xmax])

raw_input("Press Enter to exit...")
