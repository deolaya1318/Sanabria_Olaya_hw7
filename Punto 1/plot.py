import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

data = np.loadtxt(open('string_rho.dat', 'r'))
x = np.linspace(0, 100, 101)
t = np.linspace(0, 120, 121)

xu, tu = np.meshgrid(x,t)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_wireframe(xu, tu, data, rstride = 4 , cstride = 4)
ax.set_title('String position vs. time and lenght', size = 'large', color = 'r')
ax.set_xlabel('Lenght', color = 'r')
ax.set_ylabel('Time (s)', color = 'r')
ax.set_zlabel('Position of the String', color = 'r')
ax.set_xlim3d(0, 100)
ax.set_ylim3d(0, 120)
ax.set_zlim3d(-1.1, 1.1)
	
plt.savefig('cuerda_rho.pdf')
plt.close()
