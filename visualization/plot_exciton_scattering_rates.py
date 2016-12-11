# # this program reads and plots the calculated exciton scattering rates
# # in CNT conduction band.
# import numpy as np
# import matplotlib.pyplot as plt
# import sys
# import itertools

# # prepare for ploting
# figure_list = list() # this will hold the list of all the figures
# plt.ion() # enables interactive mode for ploting. No plt.show() is needed!
# colors = itertools.cycle(["r", "b", "g", "k", "c", "m", "y"]) # iterative list for cycling through colors when ploting. use color=next(colors) to get the next color on the list.
# styles = itertools.cycle(["solid", "dashed", "dash_dot"]) # iterative list for cycling through line styles when ploting. use linestyle=next(styles) to get the next color on the list.

# # set some constants
# eV = 1.6e-19 #[Jouls]
# eV_to_inverse_cm = 8065.54429
# hbar = 1.054e-34 #[Jouls.second]
# vF = 1.0e6

# #################################################################################
# fig = plt.figure()
# figure_list.append(fig)
# fig.canvas.draw()

# axes = fig.add_subplot(1,2,2)
# directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_perpendicular/"
# subdir = "trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_13_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
# file = "phonon_emission_coulomb_coupling_perpendicular.Ep_singlet_to_A2_singlet.dat"
# data = np.loadtxt(directory+subdir+file, skiprows=0)
# energy = data[0,:]/eV
# scattering_rate = data[1,:]
# axes.plot(1.0e9*scattering_rate[:],energy, linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)

# min_energy = np.amin(energy)
# max_energy = np.amax(energy)

# ymin = 1.03797
# ymax = max_energy

# axes.set_ylim([ymin,ymax])

# ################################################################################

# axes = fig.add_subplot(1,2,1)

# file  = "i_exciton.dispersion.dat"
# data = np.loadtxt(directory+subdir+file, skiprows=0)
# k_vec = data[0,:]*1e-9
# exciton_energy = (data[1:-1,:]).T
# axes.plot(k_vec, exciton_energy[:,:], linewidth=5.0, linestyle="solid", color="blue", marker="", markersize=10.0)

# # file  = "m_exciton.dispersion.dat"
# # data = np.loadtxt(directory+subdir+file, skiprows=0)
# # k_vec = data[0,:]*1e-9
# # exciton_energy = (data[1:-1,:]).T
# # axes.plot(k_vec, exciton_energy[:,:], linewidth=5.0, linestyle="solid", color="black", marker="", markersize=10.0)

# file  = "f_exciton.dispersion.dat"
# data = np.loadtxt(directory+subdir+file, skiprows=0)
# k_vec = data[0,:]*1e-9
# exciton_energy = (data[1:-1,:]).T
# axes.plot(k_vec, exciton_energy[:,:], linewidth=5.0, linestyle="solid", color="red", marker="", markersize=10.0)

# axes.set_xlim([-0.5,0.5])
# axes.set_ylim([ymin,ymax])



# # plot phonon dispersions
# # ################################################################################
# cnt_name = "cnt1"
# k_vec = np.loadtxt(directory+subdir+cnt_name+".phonon_k_vector.dat", skiprows=0)
# k_vec = k_vec*1.e-9

# # filename = directory+cnt_name+".phonon_energy.dat"
# phonon_energy = np.loadtxt(directory+subdir+cnt_name+".phonon_energy.dat", skiprows=0)
# phonon_energy = phonon_energy.T
# phonon_energy = phonon_energy/eV
# # phonon_energy = min_energy+phonon_energy
# phonon_energy = 1.03797+phonon_energy

# for i in range(0,6):
# 	axes.plot(k_vec,phonon_energy[:,i], linewidth=2.0, color="black", linestyle="solid")

# fig.savefig(directory+"figures/10Ep_to_13A2_phonon_emission_coulomb_coupling_perpendicular.pdf")

# input("Press Enter to exit...")



#####################################################################################################################################
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

#################################################################################
fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(2,1,1)
directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_parallel/"
subdir = "trf_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
file = "phonon_emission_coulomb_coupling_parallel.A2_singlet_to_A2_singlet.dat"
data = np.loadtxt(directory+subdir+file, skiprows=0)
energy = data[0,:]/eV
energy = energy-energy[0]
energy = np.exp(-energy/0.025)
energy = energy/np.sum(energy)
scattering_rate = data[1,:]
axes.plot(energy,1.0e0*scattering_rate[:], linewidth=5.0, color="red", linestyle="solid", marker="o", markersize=10.0)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_parallel/"
subdir = "trf_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_11_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
file = "phonon_emission_coulomb_coupling_parallel.A2_singlet_to_A2_singlet.dat"
data = np.loadtxt(directory+subdir+file, skiprows=0)
energy = data[0,:]/eV
energy = energy-energy[0]
energy = np.exp(-energy/0.025)
energy = energy/np.sum(energy)
scattering_rate = data[1,:]
axes.plot(energy,1.0e0*scattering_rate[:], linewidth=5.0, color="blue", linestyle="solid", marker="o", markersize=10.0)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_parallel/"
subdir = "trf_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_13_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
file = "phonon_emission_coulomb_coupling_parallel.A2_singlet_to_A2_singlet.dat"
data = np.loadtxt(directory+subdir+file, skiprows=0)
energy = data[0,:]/eV
energy = energy-energy[0]
energy = np.exp(-energy/0.025)
energy = energy/np.sum(energy)
scattering_rate = data[1,:]
axes.plot(energy,1.0e0*scattering_rate[:], linewidth=5.0, color="black", linestyle="solid", marker="o", markersize=10.0)

min_energy = np.amin(energy)
max_energy = np.amax(energy)

axes.set_xlim([min_energy,max_energy])
axes.set_ylim([1e10,1e13])
# axes.set_xscale('log')
axes.set_yscale('log')

axes = fig.add_subplot(2,1,2)
directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_parallel/"
subdir = "trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
file = "phonon_emission_coulomb_coupling_parallel.Ep_singlet_to_A2_singlet.dat"
data = np.loadtxt(directory+subdir+file, skiprows=0)
energy = data[0,:]/eV
energy = energy-energy[0]
energy = np.exp(-energy/0.025)
energy = energy/np.sum(energy)
scattering_rate = data[1,:]
axes.plot(energy,1.0e0*scattering_rate[:], linewidth=5.0, color="red", linestyle="solid", marker="o", markersize=10.0)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_parallel/"
subdir = "trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_11_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
file = "phonon_emission_coulomb_coupling_parallel.Ep_singlet_to_A2_singlet.dat"
data = np.loadtxt(directory+subdir+file, skiprows=0)
energy = data[0,:]/eV
energy = energy-energy[0]
energy = np.exp(-energy/0.025)
energy = energy/np.sum(energy)
scattering_rate = data[1,:]
axes.plot(energy,1.0e0*scattering_rate[:], linewidth=5.0, color="blue", linestyle="solid", marker="o", markersize=10.0)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/phonon_emission_coulomb_coupling_parallel/"
subdir = "trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_13_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
file = "phonon_emission_coulomb_coupling_parallel.Ep_singlet_to_A2_singlet.dat"
data = np.loadtxt(directory+subdir+file, skiprows=0)
energy = data[0,:]/eV
energy = energy-energy[0]
energy = np.exp(-energy/0.025)
energy = energy/np.sum(energy)
scattering_rate = data[1,:]
axes.plot(energy,1.0e0*scattering_rate[:], linewidth=5.0, color="black", linestyle="solid", marker="o", markersize=10.0)

min_energy = np.amin(energy)
max_energy = np.amax(energy)

axes.set_xlim([min_energy,max_energy])
axes.set_ylim([1e10,1e13])
# axes.set_xscale('log')
axes.set_yscale('log')

fig.savefig(directory+"figures/parallel_all_coulomb_coupling_phonon_emission.pdf")

input("Press Enter to exit...")