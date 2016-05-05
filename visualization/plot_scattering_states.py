# this program reads and plots the electronic and phononic states that conserve both momentum and energy in an electron-phonon scattering process.
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
# load all the needed data and calculate the needed numbers such as Nu and ikc_max

# load the list of scattering states
scattering_state_list = np.loadtxt(directory+cnt_name+".electron_phonon_scattering_states.dat", skiprows=0)

# load the electronic dispersions
k_vec = np.loadtxt(directory+cnt_name+".electron_k_vector.dat", skiprows=0)
ikc_max = (k_vec.shape[0]-1)/2

filename = directory+cnt_name+".electron_conduction_band.dat"
E_conduction = np.loadtxt(filename, skiprows=0)
E_conduction = E_conduction.T
E_conduction = E_conduction/eV

Nu = E_conduction.shape[1]

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

fig = plt.figure()
figure_list.append(fig)
axes = fig.add_subplot(111)
fig.canvas.draw()

for i in range(0,E_conduction.shape[1]):
	axes.plot(k_vec,E_conduction[:,i], linewidth=2.0, linestyle="solid")
	axes.plot(k_vec,E_valence[:,i], linewidth=2.0, linestyle="solid")

for i in range(0,phonon_energy.shape[1]):
	axes.plot(q_vec,phonon_energy[:,i], linewidth=2.0, linestyle="solid")

xmin = min(k_vec)
xmax = max(k_vec)
axes.set_xlim([xmin,xmax])

ymin = 1.01*np.amin(E_valence)
ymax = 1.01*np.amax(E_conduction)
axes.set_ylim([ymin,ymax])

################################################################################
# plot the energy and momentum conserving electronic states.

fig = plt.figure()
figure_list.append(fig)
axes = fig.add_subplot(111)
fig.canvas.draw()

for i in range(0,scattering_state_list.shape[0]):
	scattering_state_number = scattering_state_list[i,0]
	ib = scattering_state_list[i,1]
	mu_e = scattering_state_list[i,2]
	mu_e_2 = scattering_state_list[i,3]
	mu_ph = scattering_state_list[i,4]
	ik_e = scattering_state_list[i,5]
	ik_e_2 = scattering_state_list[i,6]
	iq_ph = scattering_state_list[i,7]

	print "scattering_state_number = ", scattering_state_number

	axes.clear()

	axes.plot(k_vec,E_conduction[:,int(mu_e + Nu/2-1)], linewidth=2.0, linestyle="solid", marker='o',  markersize=5)
	axes.plot(k_vec,E_conduction[:,int(mu_e_2+Nu/2-1)], linewidth=2.0, linestyle="solid", marker='o',  markersize=5)
	axes.plot(q_vec+k_vec[int(ik_e + ikc_max)],E_conduction[int(ik_e + ikc_max),int(mu_e + Nu/2-1)]-phonon_energy[:,int((ib-1)*(2*Nu-1)+(mu_ph+Nu-1))], linewidth=2.0, linestyle="solid", marker='o',  markersize=5)

	axes.plot(k_vec[int(ik_e + ikc_max)],E_conduction[int(ik_e + ikc_max),int(mu_e + Nu/2-1)], marker='o', color="red",  markersize=15)
	axes.plot(k_vec[int(ik_e_2+ikc_max)],E_conduction[int(ik_e_2+ikc_max),int(mu_e_2+Nu/2-1)], marker='o', color="blue",  markersize=15)

	xmin = min(k_vec)
	xmax = max(k_vec)
	axes.set_xlim([xmin,xmax])

	ymin = 1.01*np.amin(E_conduction)
	ymax = 1.01*np.amax(E_conduction)
	axes.set_ylim([ymin,ymax])

	raw_input("Press Enter to exit...")
