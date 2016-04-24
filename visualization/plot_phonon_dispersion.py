# this python program reads the calculated phonon dispersion of carbon nanotubes and plots the dispersion curves.
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


directory = "/home/amirhossein/research/exciton/data/transfer_rates/tmp_005/"
cnt_name = "cnt1"

################################################################################
k_vec = np.loadtxt(directory+cnt_name+".phonon_k_vector.dat", skiprows=0)

filename = directory+cnt_name+".phonon_energy.dat"
phonon_energy = np.loadtxt(filename, skiprows=0)
phonon_energy = phonon_energy.T

fig = plt.figure()
figure_list.append(fig)
axes = fig.add_subplot(111)

for i in range(0,phonon_energy.shape[1]):
	axes.plot(k_vec,phonon_energy[:,i], linewidth=2.0, linestyle="solid")

xmin = min(k_vec)
xmax = max(k_vec)
axes.set_xlim([xmin,xmax])

ymin = np.amin(phonon_energy)
ymax = np.amax(phonon_energy)
axes.set_ylim([ymin,ymax])

fig.canvas.draw()

raw_input("Press Enter to continue...")
