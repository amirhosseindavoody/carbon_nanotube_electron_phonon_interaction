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

directory = "/home/amirhossein/research/exciton/data/transfer_rates/1_transfer_10_0_Ep_singlet_iSub_1_length_0nm_center_0nm_Ckappa_2.0_to_11_0_A2_singlet_iSub_1_length_0nm_center_0nm_Ckappa_2.0_C2C_1.2nm_1.2nm_theta_0_90/"
cnt_name = "cnt1"
exciton_name = ("A1_singlet", "A2_singlet", "Ep_singlet", "Em_singlet")

################################################################################
# load the calculated scattering rates for emission and absorption process

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

xmin = +100
xmax = -100

ymin = 1e0
ymax = 1e7

axes = fig.add_subplot(1,1,1)

for i in range(0,3):
	# axes = fig.add_subplot(2,2,i+1)
	#
	# xmin = +100
	# xmax = -100

	for j in range(1,2):

		data = np.loadtxt(directory+"phonon_assisted_scattering_rate_emission."+exciton_name[i]+"_to_"+exciton_name[j]+".dat", skiprows=0)
		energy_mesh_emission = data[0,:]/eV
		# energy_mesh_emission = energy_mesh_emission - np.amin(energy_mesh_emission)
		hopping_rate_emission = data[1,:].T
		axes.plot(energy_mesh_emission, hopping_rate_emission[:], linewidth=5.0, linestyle="solid", marker="", markersize=10.0)

		xmin = min(xmin,np.amin(energy_mesh_emission))
		xmax = max(xmax,np.amax(energy_mesh_emission))
		axes.set_xlim([xmin,xmax])

		axes.set_yscale('log')
		# axes.set_ylim([ymin,ymax])

		input(exciton_name[i]+" to "+exciton_name[j])

input("Press Enter to exit...")
