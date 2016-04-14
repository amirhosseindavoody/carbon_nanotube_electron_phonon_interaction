import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.pyplot import plot, ion, show
# ion() # enables interactive mode
# plot([1,2,3]) # result shows immediatelly (implicit draw())
#
# print 'continue computation'
# raw_input("Press Enter to continue...")
#
# # at the end call show to ensure window won't close.
# show()
#
# raw_input("Press Enter to continue...")


x = np.arange(0.0, 10.0, 0.1)
plt.ion()
plt.plot(x,x)
# plt.draw()
print 'continue computation'

raw_input("Press Enter to continue...")
# plt.show()

# y = pow(x,2)
y=x**2
plt.plot(x,y)
# plt.draw()
# print 'continue computation'

# raw_input("Press Enter to continue...")

# at the end call show to ensure window won't close.
# plt.show()
raw_input("Press Enter to continue...")
