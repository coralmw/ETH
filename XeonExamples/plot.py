import numpy as np
import matplotlib.pyplot as plt

upup = np.loadtxt("results/upup")
downdown = np.loadtxt("results/downdown")


plt.plot(upup)
plt.plot(downdown)
plt.plot(downdown+upup)


plt.show()
