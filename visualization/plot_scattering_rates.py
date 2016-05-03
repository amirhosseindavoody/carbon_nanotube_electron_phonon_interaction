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

directory = "/home/amirhossein/research/exciton/data/transfer_rates/tmp_001/"
cnt_name = "cnt1"

################################################################################
# load the calculated scattering rates for emission and absorption process
data = np.loadtxt(directory+cnt_name+".electron_phonon_scattering_rate_emission.dat", skiprows=0)
energy_mesh_emission = data[0,:]/eV
scattering_rate_emission = data[1,:]

data = np.loadtxt(directory+cnt_name+".electron_phonon_scattering_rate_absorption.dat", skiprows=0)
energy_mesh_absorption = data[0,:]/eV
scattering_rate_absorption = data[1,:]

energy_mesh_emission = energy_mesh_emission - np.amin(energy_mesh_emission)
energy_mesh_absorption = energy_mesh_absorption - np.amin(energy_mesh_absorption)

################################################################################
# plot the electron-phonon scattering rates
fig = plt.figure()
figure_list.append(fig)
axes = fig.add_subplot(1,2,1)
fig.canvas.draw()

axes.plot(scattering_rate_emission, energy_mesh_emission, color="blue", linewidth=3.0, linestyle="solid", marker="")
axes.plot(scattering_rate_absorption, energy_mesh_absorption, color="red", linewidth=3.0, linestyle="solid", marker="")

# xmin = 1.e11
# xmax = 1.e15
# axes.set_xlim([xmin,xmax])

ymin = min(energy_mesh_emission)
ymax = max(energy_mesh_emission)
axes.set_ylim([ymin,ymax])

axes.set_xscale('log')

################################################################################
# load the electronic dispersions
k_vec = np.loadtxt(directory+cnt_name+".electron_k_vector.dat", skiprows=0)

filename = directory+cnt_name+".electron_conduction_band.dat"
E_conduction = np.loadtxt(filename, skiprows=0)
E_conduction = E_conduction.T
E_conduction = E_conduction/eV

E_conduction = E_conduction - np.amin(E_conduction)

filename = directory+cnt_name+".electron_valence_band.dat"
E_valence = np.loadtxt(filename, skiprows=0)
E_valence = E_valence.T
E_valence = E_valence/eV

# load the phononic dispersion
q_vec = np.loadtxt(directory+cnt_name+".phonon_k_vector.dat", skiprows=0)

filename = directory+cnt_name+".phonon_energy.dat"
phonon_energy = np.loadtxt(filename, skiprows=0)
phonon_energy = phonon_energy.T
phonon_energy = phonon_energy/eV

################################################################################
# plot the electroic dispersion and phononic dispersions all in one figure.

axes = fig.add_subplot(1,2,2)
fig.canvas.draw()

for i in range(0,E_conduction.shape[1]):
	axes.plot(k_vec,E_conduction[:,i], linewidth=3.0, linestyle="solid")

for i in range(0,phonon_energy.shape[1]):
	axes.plot(q_vec,np.amin(E_conduction)+phonon_energy[:,i], linewidth=2.0, linestyle="solid")

xmin = min(k_vec)
xmax = max(k_vec)
axes.set_xlim([xmin,xmax])

ymin = min(energy_mesh_emission)
ymax = max(energy_mesh_emission)
axes.set_ylim([ymin,ymax])


raw_input("Press Enter to exit...")
