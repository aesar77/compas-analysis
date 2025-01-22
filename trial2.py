import os, sys
import numpy as np               # for handling arrays
import numpy.ma as ma
import pandas as pd
import h5py as h5                # for reading the COMPAS data
import time                      # for finding computation time
from datetime import datetime
import matplotlib.pyplot as plt  #for plotting
import matplotlib

i = np.random.uniform(0,180,100000000)
cosiDC2 = (np.abs(np.cos(i)) <= 2.22819667e-6)
cosiDC2arc = (i <= 2*(90 - 89.99987233))

print(i)
print(cosiDC2)
print(len(i))
print(sum(cosiDC2))
print(sum(cosiDC2arc))

# plt.scatter(i, np.cos(i))

# plt.show()
# i = i*cosiDC2

# plt.scatter(i, np.cos(i))

# plt.show()