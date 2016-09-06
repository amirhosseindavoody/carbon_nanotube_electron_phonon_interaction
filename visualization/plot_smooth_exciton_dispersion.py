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

axes = fig.add_subplot(1,1,1)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/trf_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_11_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"

filename = directory+"cnt2.exciton_dispersion.dat"
ex_energy = np.loadtxt(filename, skiprows=0)
ex_energy = np.transpose(ex_energy)

# build the k_vec with the appropriate size of nKcm
nKcm = ex_energy.shape[0]
k_vec_exciton = np.linspace(-(nKcm-1)/2,+(nKcm-1)/2,num=nKcm,endpoint=True)

nX = ex_energy.shape[1]
print("nKcm = ", nKcm)
print("nX = ", nX)

# arr = ex_energy[:,1]
# arr = np.flip(arr)
# print(arr.shape)

# input("enter")

# for iX in range(0,nX):
axes.plot(k_vec_exciton, ex_energy[:,0], linewidth=3.0, linestyle="solid", color="red", marker="")
axes.plot(k_vec_exciton, ex_energy[:,1], linewidth=3.0, linestyle="solid", color="blue", marker="")
# ex_energy = np.flipud(ex_energy)
# axes.plot(k_vec_exciton, ex_energy[:,1], linewidth=3.0, linestyle="solid", color="red", marker="")


# xmin = np.amin(k_vec_exciton)
# xmax = np.amax(k_vec_exciton)
xmin=-70
xmax = -xmin
axes.set_xlim([xmin,xmax])

ymin = np.amin(ex_energy)
ymax = ymin + 0.18
axes.set_ylim([ymin,ymax])

input("Press Enter to exit...")
