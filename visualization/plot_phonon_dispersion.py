# this program reads and plots the calculated phonon dispersion of carbon nanotubes.
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
# load k-vector and phonon energy dispersion data.
k_vec = np.loadtxt(directory+cnt_name+".phonon_k_vector.dat", skiprows=0)

filename = directory+cnt_name+".phonon_energy.dat"
phonon_energy = np.loadtxt(filename, skiprows=0)
phonon_energy = phonon_energy.T
phonon_energy = phonon_energy/eV

################################################################################
# plot all phonon dispersion energies
fig = plt.figure()
figure_list.append(fig)
axes = fig.add_subplot(111)
fig.canvas.draw()

for i in range(0,phonon_energy.shape[1]):
	axes.plot(k_vec,phonon_energy[:,i], linewidth=2.0, linestyle="solid")

xmin = min(k_vec)
xmax = max(k_vec)
axes.set_xlim([xmin,xmax])

ymin = 1.1*np.amin(phonon_energy)
ymax = 1.1*np.amax(phonon_energy)
axes.set_ylim([ymin,ymax])

# raw_input("Press Enter to continue...")

################################################################################
# plot phonon dispersion for a specific phonon branch (ib) and mu_ph

Nu = (phonon_energy.shape[1]/6+1)/2

print "ib should be in interval [ 1 , 6 ]"
ib = int(raw_input("Input phonon branch index: "))
print "mu_ph should be in interval [", 1-Nu, ",", Nu-1,"]"
mu_ph = int(raw_input("Input mu_ph value: "))

print "ib = ", ib
print "mu_ph = ", mu_ph

i = (ib-1)*(2*Nu-1)+(mu_ph+Nu-1)
axes.plot(k_vec,phonon_energy[:,i], linewidth=10, linestyle="solid", color="red")

raw_input("Press Enter to exit...")
