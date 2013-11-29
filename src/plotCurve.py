import numpy as np
#import matplotlib.pyplot as plt
import pylab as pl
import sys

x,y = np.loadtxt(sys.argv[1], unpack = True)

#xx, xy = np.loadtxt(sys.argv[2], unpack = True)

pl.clf()
pl.plot(x[:51], y[:51]/50, 'b^')
#plt.plot(xx, xy, 'r')

pl.xlabel('n')
pl.ylabel(r'$f_n$')
#pl.title('f(x) vs x')
#pl.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
#plt.grid(True)
pl.savefig('complexFFT.pdf')


