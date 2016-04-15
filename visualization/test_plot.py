import matplotlib.pyplot as plt
import numpy as np

plt.ion()

figure_list = list()

x = np.arange(0,5,0.01)
y = x**2/25.0
fig = plt.figure()
figure_list.append(fig)
ax = figure_list[0].add_subplot(111)
ax.plot(x, y)
# plt.plot(x,y)
# ax1 = plt.gca()

z = np.sin(x)
fig = plt.figure()
figure_list.append(fig)
ax = figure_list[1].add_subplot(111)
ax.plot(x, z)
# plt.plot(x,z)
# ax2 = plt.gca()

w = np.cos(x)
ax = figure_list[0].add_subplot(111)
ax.plot(x, w) # can continue plotting on the first axis
figure_list[0].canvas.draw()
# ax1.plot(x,w)

ymin = -1
ymax = 2
ax.set_ylim([ymin,ymax])

# plt.show()

raw_input("press enter to continue...")
