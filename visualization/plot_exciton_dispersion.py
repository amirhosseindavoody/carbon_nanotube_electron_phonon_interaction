# this program reads and plots the calculated exciton energy dispersions
# for a carbon nanotube. The electron dispersion is also ploted for comparision.
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

################################################################################
# load the exciton dispersions
directory = "/home/amirhossein/research/exciton/data/exciton_dispersion/CNT(10,00)-nkg(1001)-nr(0200)-E_th(0.5)-Kcm_max(1.5)-i_sub(1)-Ckappa(1.0)/"
k_vec_temporary = np.loadtxt(directory+"kVec_fine.dat", skiprows=0)
dk = k_vec_temporary[1]-k_vec_temporary[0]
del k_vec_temporary

filename = directory+"Ex0_A2.dat"
Ex0_A2 = np.loadtxt(filename, skiprows=0)
Ex0_A2 = Ex0_A2/eV
Ex0_A2 = Ex0_A2-np.amin(Ex0_A2)

# build the k_vec with the appropriate size of nKcm
nKcm = Ex0_A2.shape[0]
k_vec_exciton = np.linspace(-(nKcm-1)/2*dk,+(nKcm-1)/2*dk,num=nKcm,endpoint=True)

################################################################################
# load the electron dispersions for conduction band
directory = "/home/amirhossein/research/exciton/data/transfer_rates/final_result/"
cnt_name = "cnt1"
k_vec_electron = np.loadtxt(directory+cnt_name+".electron_k_vector.dat", skiprows=0)

filename = directory+cnt_name+".electron_conduction_band.dat"
Ec = np.loadtxt(filename, skiprows=0)
Ec = Ec.T
Ec = Ec/eV
Ec = Ec-np.amin(Ec)

################################################################################
# plot the excitonic dispersion.

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(1,1,1)

for i in range(0,Ex0_A2.shape[1]):
	axes.plot(k_vec_exciton, Ex0_A2[:,i], linewidth=3.0, linestyle="solid")


xmin = min(k_vec_exciton)
xmax = max(k_vec_exciton)
axes.set_xlim([xmin,xmax])

ymin = np.amin(Ex0_A2)
ymax = np.amax(Ex0_A2)
axes.set_ylim([ymin,ymax])

################################################################################
# plot the electron dispersions

for i in range(0,Ec.shape[1]):
	axes.plot(k_vec_electron,Ec[:,i], linewidth=3.0, linestyle="solid", color="black")


raw_input("Press Enter to exit...")
