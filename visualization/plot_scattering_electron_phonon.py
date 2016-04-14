# this python program reads the calculated phonon dispersion of carbon nanotubes and plots the dispersion curves.
import numpy as np
import matplotlib.pyplot as plt
import sys
import itertools

eV = 1.6e-19 #[Jouls]
eV_to_inverse_cm = 8065.54429
hbar = 1.054e-34 #[Jouls.second]
vF = 1.0e6

directory = "/home/amirhossein/research/exciton/data/transfer_rates/tmp_002/"

################################################################################
k_vec = np.loadtxt(directory+"cnt1.electron_k_vector.dat", skiprows=0)

colors = itertools.cycle(["r", "b", "g", "k", "c", "m", "y"])
# styles = itertools.cycle(["solid", "dashed", "dash_dot", "dotted"])
styles = itertools.cycle(["solid", "dashed", "dash_dot"])

for file_number in range(1,6):
	filename = directory+"cnt1.electron_phonon_matrix_element_branch_"+str(file_number)+".dat"

	matrix_element_branch = np.loadtxt(filename, skiprows=0)
	matrix_element_branch = matrix_element_branch.T

	line_style = next(styles)

	if (matrix_element_branch.ndim > 1):
		for i in range(0,matrix_element_branch.shape[1]):
			plt.plot(k_vec,matrix_element_branch[:,i], color=next(colors), linewidth=2.0, linestyle="solid")

		axes = plt.gca()
		xmin = min(k_vec)
		xmax = max(k_vec)
		ymin = -1.0e-20
		ymax = 8.0e-19
		axes.set_xlim([xmin,xmax])
		axes.set_ylim([ymin,ymax])
		plt.show()
