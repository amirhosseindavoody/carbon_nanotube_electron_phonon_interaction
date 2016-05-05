# this program reads and plots the calculated electron scattering rates in CNT conduction band.
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
# load the electronic and phononic dispersions
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
# load the calculated scattering rates for emission and absorption process
data = np.loadtxt(directory+cnt_name+".electron_phonon_scattering_rate_emission.dat", skiprows=0)
energy_mesh_emission = data[0,:]/eV
scattering_rate_emission = data[1:,:]
scattering_rate_emission = scattering_rate_emission.T

data = np.loadtxt(directory+cnt_name+".electron_phonon_scattering_rate_absorption.dat", skiprows=0)
energy_mesh_absorption = data[0,:]/eV
scattering_rate_absorption = data[1:,:]
scattering_rate_absorption = scattering_rate_absorption.T

energy_mesh_emission = energy_mesh_emission - np.amin(energy_mesh_emission)
energy_mesh_absorption = energy_mesh_absorption - np.amin(energy_mesh_absorption)

################################################################################
# plot the electroic dispersion and phononic dispersions all in one figure.

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(3,1,1)

for i in range(0,E_conduction.shape[1]):
	axes.plot(E_conduction[:,i], k_vec, linewidth=3.0, linestyle="solid")

for i in range(0,phonon_energy.shape[1]):
	axes.plot(np.amin(E_conduction)+phonon_energy[:,i], q_vec, linewidth=2.0, linestyle="solid")

xmin = min(energy_mesh_emission)
xmax = max(energy_mesh_emission)
axes.set_xlim([xmin,xmax])

ymin = min(k_vec)
ymax = max(k_vec)
axes.set_ylim([ymin,ymax])

################################################################################
# plot the electron-phonon scattering rates

axes_1 = fig.add_subplot(3,1,2)
axes_2 = fig.add_subplot(3,1,3)

scattering_rate_emission_total = np.zeros(scattering_rate_emission.shape[0])
scattering_rate_absorption_total = np.zeros(scattering_rate_absorption.shape[0])

for i in range(0,scattering_rate_emission.shape[1]):
	scattering_rate_emission_total += scattering_rate_emission[:,i]
	scattering_rate_absorption_total += scattering_rate_absorption[:,i]

axes_1.plot(energy_mesh_emission, scattering_rate_emission_total, color="black", linewidth=5.0, linestyle="solid", marker="")
axes_2.plot(energy_mesh_absorption, scattering_rate_absorption_total, color="black", linewidth=5.0, linestyle="solid", marker="")

for i in range(0,scattering_rate_emission.shape[1]):
	axes_1.plot(energy_mesh_emission, scattering_rate_emission[:,i], linewidth=3.0, linestyle="solid", marker="")
	axes_2.plot(energy_mesh_absorption, scattering_rate_absorption[:,i], linewidth=3.0, linestyle="solid", marker="")

	xmin = min(energy_mesh_emission)
	xmax = max(energy_mesh_emission)
	axes_1.set_xlim([xmin,xmax])

	axes_1.set_yscale('log')

	xmin = min(energy_mesh_emission)
	xmax = max(energy_mesh_emission)
	axes_2.set_xlim([xmin,xmax])

	axes_2.set_yscale('log')

	scattering_rate_emission_total += scattering_rate_emission[:,i]
	scattering_rate_absorption_total += scattering_rate_absorption[:,i]

	print "ib = ", i
	raw_input("Press Enter to continue...")


raw_input("Press Enter to exit...")
