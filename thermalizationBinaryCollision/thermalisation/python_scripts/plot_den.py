import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys

filename = 'results_1000.txt'
path = './data/files/'

NUM_TS = 12000
write_interval = 100
NC = 400
ni = NC+1


DATA_TS = int(NUM_TS/write_interval)+1

if os.path.exists(pjoin(path,filename)):
    x, ndA, ndB = np.loadtxt(pjoin(path, filename), unpack=True)
else:
    print('No Data')
    exit()

x = x.reshape(DATA_TS, (NC+1))
ndA = ndA.reshape(DATA_TS, (NC+1))
ndB = ndB.reshape(DATA_TS, (NC+1))

plt.figure(1)
for i in range(DATA_TS): #DATA_TS
    plt.cla()
    plt.plot(x[i,:], ndA[i,:])
    plt.plot(x[i,:], ndB[i,:])
    plt.xlabel('x')
    plt.ylabel('neutral density')
    plt.ylim(0,4E17)
    plt.title(i)
    plt.savefig("./figure_den/den_%d.png"%(i))
    plt.pause(0.05)

#plt.savefig('plot_den.png')
plt.show()
