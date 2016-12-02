import numpy as np
import matplotlib.pyplot as plt

res = open("results",'r')

rhoreduced = np.zeros((4, 4, 1000))

for line in res:
  idx1, idx2, *data = line.split(' ')
  rhoreduced[int(idx1),int(idx2)] = np.array(data[:-1])

f, axarr = plt.subplots(2, sharex=True)
for traceidx in range(4):
    axarr[0].plot(rhoreduced[traceidx, traceidx])

axarr[0].plot(np.trace(rhoreduced)) # total probability

axarr[1].plot(np.sum(rhoreduced, axis=(0,1)) - np.trace(rhoreduced))

plt.show()
