import numpy as np
import matplotlib.pyplot as plt
import sys

with open(sys.argv[1]) as f: 
    val = [float(line) for line in f]

plt.hist(val,50)

plt.show()
