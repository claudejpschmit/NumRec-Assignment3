import numpy as np
import matplotlib.pyplot as plt
import sys

x1,y1 = np.loadtxt(sys.argv[1], unpack = True)
plt.plot(x1,y1)
plt.xlabel(r'$x$')
plt.ylabel(r'$F(x)$')
plt.title(r'$F(x) vs x$')
plt.grid(True)
plt.ylim([-1,6])
plt.show()
