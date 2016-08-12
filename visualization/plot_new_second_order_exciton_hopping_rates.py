# this program reads and plots the calculated exciton energy dispersions
# for a carbon nanotube. The electron dispersion is also ploted for comparision.
import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools

# prepare for ploting
figure_list = list() # this will hold the list of all the figures
plt.ion() # enables interactive mode for ploting. No plt.show() is needed!

# set some constants
eV = 1.6e-19 #[Jouls]
eV_to_inverse_cm = 8065.54429
hbar = 1.054e-34 #[Jouls.second]
vF = 1.0e6

################################################################################
# plot the excitonic dispersion.

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(1,2,1)

directory = "/home/amirhossein/research/exciton/data/exciton_dispersion/CNT(17,00)-nkg(1001)-nr(0200)-E_th(0.5)-Kcm_max(1.5)-i_sub(1)-Ckappa(2.0)/"
k_vec_temporary = np.loadtxt(directory+"kVec_fine.dat", skiprows=0)
dk = k_vec_temporary[1]-k_vec_temporary[0]
del k_vec_temporary

filename = directory+"Ex0_A2.dat"
Ex0_A2 = np.loadtxt(filename, skiprows=0)
Ex0_A2 = Ex0_A2/eV

# build the k_vec with the appropriate size of nKcm
nKcm = Ex0_A2.shape[0]
k_vec_exciton = np.linspace(-(nKcm-1)/2*dk,+(nKcm-1)/2*dk,num=nKcm,endpoint=True)

for i in range(0,Ex0_A2.shape[1]):
	axes.plot(k_vec_exciton, Ex0_A2[:,i], linewidth=3.0, linestyle="solid", color="red", marker="")

del Ex0_A2
del k_vec_exciton


directory = "/home/amirhossein/research/exciton/data/exciton_dispersion/CNT(10,00)-nkg(1001)-nr(0200)-E_th(0.5)-Kcm_max(1.5)-i_sub(1)-Ckappa(2.0)/"
k_vec_temporary = np.loadtxt(directory+"kVec_fine.dat", skiprows=0)
dk = k_vec_temporary[1]-k_vec_temporary[0]
del k_vec_temporary

filename = directory+"Ex0_A2.dat"
Ex0_A2 = np.loadtxt(filename, skiprows=0)
Ex0_A2 = Ex0_A2/eV

# build the k_vec with the appropriate size of nKcm
nKcm = Ex0_A2.shape[0]
k_vec_exciton = np.linspace(-(nKcm-1)/2*dk,+(nKcm-1)/2*dk,num=nKcm,endpoint=True)

for i in range(0,Ex0_A2.shape[1]):
	axes.plot(k_vec_exciton, Ex0_A2[:,i], linewidth=3.0, linestyle="solid", color="blue", marker="")

filename = directory+"Ex0_A2.dat"
Ex0_A2 = np.loadtxt(filename, skiprows=0)
Ex0_A2 = Ex0_A2/eV

for i in range(0,Ex0_A2.shape[1]):
	axes.plot(k_vec_exciton, Ex0_A2[:,i], linewidth=3.0, linestyle="solid", color="black", marker="")

# xmin = min(k_vec_exciton)
# xmax = max(k_vec_exciton)
xmin = -0.5e9
xmax = +0.5e9
axes.set_xlim([xmin,xmax])

ymin = np.amin(Ex0_A2)
ymax = np.amax(Ex0_A2)
# axes.set_ylim([ymin,ymax])
axes.set_ylim([ymin,ymin+0.18])


################################################################################
# plot the second-order exciton transfer rate

axes = fig.add_subplot(1,2,2)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/"
data_temporary = np.loadtxt(directory+"Ep_to_A2.dat", skiprows=0)
energy = data_temporary[0,:]
transfer_rate = data_temporary[1,:]
del data_temporary

energy = energy + ymin

axes.plot(transfer_rate, energy, linewidth=3.0, linestyle="solid", color="blue", marker="")
axes.set_xscale('log')

# axes.set_xlim([1e1,1e7])
axes.set_ylim([ymin,ymin+0.18])

input("Press Enter to exit...")
