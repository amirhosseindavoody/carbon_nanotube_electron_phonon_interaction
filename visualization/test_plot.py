from matplotlib.pyplot import plot, ion, show
ion() # enables interactive mode
plot([1,2,3]) # result shows immediatelly (implicit draw())

print 'continue computation'

# at the end call show to ensure window won't close.
show()
