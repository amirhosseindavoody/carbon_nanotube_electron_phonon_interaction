# this python program reads the calculated phonon dispersion of carbon nanotubes and plots the dispersion curves.
import numpy as np
import matplotlib.pyplot as plt

eV = 1.6e-19
eV_to_inverse_cm = 8065.54429

data = np.loadtxt("/home/amirhossein/research/exciton/data/transfer_rates/tmp_001/phonon_dispersion.dat", skiprows=0)
# data = np.loadtxt("/home/amirhossein/research/exciton/data/electron_phonon/CNT(10,10)-nkg(1001)-nr(0200)-E_th(0.5)-Kcm_max(1.5)-i_sub(2)-Ckappa(1.0)/phonon_dispersion.dat", skiprows=0)

k_vec = np.array(data[0,:])
omega = np.array(np.transpose(data[1:-1,:]))

# convert the units to electron-volt
omega = omega/eV

# convert the units to inverse cm.
omega = omega*eV_to_inverse_cm

plt.plot(k_vec[:],omega[:,:],'b-')
# plt.tight_layout()
plt.show()

# k_vec = data[0,:]
# omega = data[]
# plt.plot(data_dw[:,0], data_dw[:,1], 'r-')
#
# data_up = np.loadtxt("ld1.wfc.up", skiprows=1)
# plt.plot(data_up[:,0], data_up[:,1], 'b--')
#
# xmin = 0
# xmax = 10
# ymin = -0.7
# ymax = 0.3
# axes = plt.gca()
# axes.set_xlim([xmin,xmax])
# axes.set_ylim([ymin,ymax])
#
# plt.show()

# array = np.zeros((2,3),np.float)
# print "array = \n", array
#
# print "range = ", range(0,2)
#
# k=0
# for i in range(0,2):
# 	for j in range(0,3):
# 		k = k+1
# 		array[i,j] = k
#
# print "array = \n", array
#
# for row in range(0,shp[0]-2):
# 	print row
# 	plt.plot(k_vec[:],omega[row,:],'r-')
#
# plt.show()
