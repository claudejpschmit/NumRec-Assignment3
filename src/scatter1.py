import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import sys

x,y,z = np.loadtxt(sys.argv[1], unpack = True)

pl.plot(x,y, 'ro')

pl.show()

