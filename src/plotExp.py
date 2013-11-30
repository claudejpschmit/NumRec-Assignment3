import numpy as np
import matplotlib.pyplot as plot
#data = np.loadtxt("plot_me.dat")
#x,y=np.loadtxt("plot_me.dat",unpack=True) #thanks warren!
#x,y =  zip(*data)
#plot.plot(x, y, linewidth=2.0)
plot.plot(*np.loadtxt("exact1.txt",unpack=True), label = "exact")
plot.plot(*np.loadtxt("euler1.txt",unpack=True), label = "Euler")
plot.plot(*np.loadtxt("mkm1.txt",unpack=True), label = "MKM")
plot.plot(*np.loadtxt("fork1.txt",unpack=True), label = "FORK")
plot.legend()
plot.show()
