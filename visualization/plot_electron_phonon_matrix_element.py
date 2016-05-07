# this program reads and plots the calculated electron-phonon interaction matrix element.
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
k_vec = np.loadtxt(directory+"cnt1.electron_k_vector.dat", skiprows=0)

for file_number in range(1,7):
	filename = directory+cnt_name+".electron_phonon_matrix_element_branch_"+str(file_number)+".dat"

	matrix_element_branch = np.loadtxt(filename, skiprows=0)
	matrix_element_branch = matrix_element_branch.T

	fig = plt.figure()
	figure_list.append(fig)
	axes = fig.add_subplot(111)

	if (matrix_element_branch.ndim > 1):
		for i in range(0,matrix_element_branch.shape[1]):
			axes.plot(k_vec,matrix_element_branch[:,i], linewidth=2.0, linestyle="solid")

		xmin = min(k_vec)
		xmax = max(k_vec)
		ymin = -1.0
		ymax = 18.0
		axes.set_xlim([xmin,xmax])
		axes.set_ylim([ymin,ymax])
		fig.canvas.draw()


raw_input("Press Enter to continue...")
