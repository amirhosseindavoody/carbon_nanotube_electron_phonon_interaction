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

# directory = "/home/amirhossein/research/exciton/data/transfer_rates/trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_11_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"

################################################################################
# load the calculated scattering rates for emission and absorption process

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(1,1,1)

# directory = "/home/amirhossein/research/exciton/data/transfer_rates/trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
# data = np.loadtxt(directory+"phonon_assisted_scattering_rate_emission.Ep_singlet_to_A2_singlet.dat", skiprows=0)
# energy_mesh_emission = data[0,:]/eV
# scattering_rate_emission = data[1,:].T
# axes.plot(energy_mesh_emission, 1.0e9*scattering_rate_emission[:], linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)

# input("10_Ep to 10_A2")

directory = "/home/amirhossein/research/exciton/data/transfer_rates/trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
data = np.loadtxt(directory+"phonon_assisted_scattering_rate_emission.Ep_singlet_to_A2_singlet.dat", skiprows=0)
energy_mesh_emission = data[0,:]/eV
scattering_rate_emission = data[1,:].T
# axes.plot(energy_mesh_emission, 1.0e9*scattering_rate_emission[:], linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)
axes.plot(1.0e9*scattering_rate_emission[:], energy_mesh_emission, linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)

min_energy = np.amin(energy_mesh_emission)
max_energy = np.amax(energy_mesh_emission)

axes.set_ylim([min_energy,max_energy])

# input("10_Ep to 11_A2")



# directory = "/home/amirhossein/research/exciton/data/transfer_rates/trf_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_11_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
# data = np.loadtxt(directory+"phonon_assisted_scattering_rate_emission.A2_singlet_to_A2_singlet.dat", skiprows=0)
# energy_mesh_emission = data[0,:]/eV
# scattering_rate_emission = data[1,:].T
# axes.plot(energy_mesh_emission, 1.0e9*scattering_rate_emission[:], linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)

# input("10_A2 to 11_A2")

# directory = "/home/amirhossein/research/exciton/data/transfer_rates/trf_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
# data = np.loadtxt(directory+"phonon_assisted_scattering_rate_emission.A2_singlet_to_A2_singlet.dat", skiprows=0)
# energy_mesh_emission = data[0,:]/eV
# scattering_rate_emission = data[1,:].T
# axes.plot(energy_mesh_emission, 1.0e9*scattering_rate_emission[:], linewidth=5.0, linestyle="solid", marker="o", markersize=10.0)

# input("10_A2 to 10_A2")

# input("Press Enter to exit...")



################################################################################
# load the calculated scattering rates for emission and absorption process

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(1,1,1)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/rr.trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
data = np.loadtxt(directory+"i_exciton.dispersion.dat", skiprows=0)
k_vec = data[0,:]*1e-9
i_exciton_energy = (data[1:-1,:])
axes.plot(k_vec, i_exciton_energy[0,:], linewidth=5.0, linestyle="solid", color="blue", marker="", markersize=10.0)

data = np.loadtxt(directory+"m_exciton.dispersion.dat", skiprows=0)
k_vec = data[0,:]*1e-9
m_exciton_energy = (data[1:-1,:])
axes.plot(k_vec, m_exciton_energy[0,:], linewidth=5.0, linestyle="solid", color="red", marker="", markersize=10.0)

data = np.loadtxt(directory+"f_exciton.dispersion.dat", skiprows=0)
k_vec = data[0,:]*1e-9
f_exciton_energy = (data[1:-1,:])
axes.plot(k_vec, f_exciton_energy[0,:], linewidth=5.0, linestyle="solid", color="green", marker="", markersize=10.0)


# plot phonon dispersions

################################################################################
# load k-vector and phonon energy dispersion data.
cnt_name = "cnt1"
k_vec = np.loadtxt(directory+cnt_name+".phonon_k_vector.dat", skiprows=0)

filename = directory+cnt_name+".phonon_energy.dat"
phonon_energy = np.loadtxt(filename, skiprows=0)
phonon_energy = phonon_energy.T
phonon_energy = phonon_energy/eV
phonon_energy = phonon_energy + min_energy

################################################################################
# plot all phonon dispersion energies
# fig = plt.figure()
# figure_list.append(fig)
# axes = fig.add_subplot(111)
# fig.canvas.draw()

for i in range(0,phonon_energy.shape[1]):
	axes.plot(k_vec,phonon_energy[:,i], linewidth=2.0, color="black", linestyle="solid")

# xmin = min(k_vec)
# xmax = max(k_vec)
# axes.set_xlim([xmin,xmax])


axes.set_xlim([-0.5,0.5])
axes.set_ylim([min_energy,max_energy])


################################################################################
# load the calculated scattering rates for emission and absorption process

fig = plt.figure()
figure_list.append(fig)
fig.canvas.draw()

axes = fig.add_subplot(1,1,1)

directory = "/home/amirhossein/research/exciton/data/transfer_rates/r.trf_10_0_Ep_singlet_sub_1_len_0_ctr_0_dk_ratio_10_to_10_0_A2_singlet_sub_1_len_0_ctr_0_dk_ratio_10_C2C_1.2_1.2_theta_0_90/"
data = np.loadtxt(directory+"i_exciton.dispersion.dat", skiprows=0)
k_vec = data[0,:]*1e-9
i_exciton_energy = (data[1:-1,:])
axes.plot(k_vec, i_exciton_energy[0,:], linewidth=5.0, linestyle="solid", color="blue", marker="", markersize=10.0)

data = np.loadtxt(directory+"m_exciton.dispersion.dat", skiprows=0)
k_vec = data[0,:]*1e-9
m_exciton_energy = (data[1:-1,:])
axes.plot(k_vec, m_exciton_energy[0,:], linewidth=5.0, linestyle="solid", color="red", marker="", markersize=10.0)

data = np.loadtxt(directory+"f_exciton.dispersion.dat", skiprows=0)
k_vec = data[0,:]*1e-9
f_exciton_energy = (data[1:-1,:])
axes.plot(k_vec, f_exciton_energy[0,:], linewidth=5.0, linestyle="solid", color="green", marker="", markersize=10.0)

axes.set_xlim([-0.5,0.5])
axes.set_ylim([min_energy,max_energy])

input("Press Enter to exit...")