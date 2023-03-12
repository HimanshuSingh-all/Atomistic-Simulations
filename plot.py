import numpy as np
import matplotlib.pyplot as plt

points = np.loadtxt('macros.txt', skiprows = 1)

plt.plot(points[:,0], points[:,1], label = 'KE')
plt.plot(points[:,0], points[:,2], label = 'U')
plt.plot(points[:,0], points[:,1]+points[:,2] , label = 'E')

plt.legend()
plt.show()
