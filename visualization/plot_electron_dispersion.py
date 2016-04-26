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

k_vec = np.loadtxt(directory+cnt_name+".electron_k_vector.dat", skiprows=0)

data = np.loadtxt("/home/amirhossein/research/exciton/data/transfer_rates/tmp_001/cnt1.electron_conduction_band.dat", skiprows=0)
Ec = np.array(np.transpose(data[:,:]))

filename = directory+cnt_name+".electron_conduction_band.dat"
Ec = np.loadtxt(filename, skiprows=0)
Ec = Ec.T
Ec = Ec/eV

filename = directory+cnt_name+".electron_valence_band.dat"
Ev = np.loadtxt(filename, skiprows=0)
Ev = Ev.T
Ev = Ev/eV

fig = plt.figure()
figure_list.append(fig)
axes = fig.add_subplot(111)

for i in range(0,Ec.shape[1]):
	axes.plot(k_vec,Ec[:,i], linewidth=2.0, linestyle="solid")
	axes.plot(k_vec,Ev[:,i], linewidth=2.0, linestyle="solid")

xmin = min(k_vec)
xmax = max(k_vec)
axes.set_xlim([xmin,xmax])

ymin = 1.01*np.amin(Ev)
ymax = 1.01*np.amax(Ec)
axes.set_ylim([ymin,ymax])

fig.canvas.draw()

raw_input("Press Enter to exit...")
