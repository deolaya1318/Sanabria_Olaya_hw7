from pylab import *
import sys

nombre = sys.argv[1];
arr = nombre.split("_")
arr1 = arr[1].split(".dat")
tiempo = arr1[0]

datos = loadtxt(nombre)

x = datos[:,0]
v = datos[:,1]
p = datos[:,2]
d = datos[:,3]

figura = figure(figsize=(10, 14))

ax = figura.add_subplot(3, 1, 1)
ax.plot(x, v)
ax.set_xlabel("$x (m)$")
ax.set_ylabel("$Vel (m/s)$")

ax = figura.add_subplot(3, 1, 2)
ax.plot(x, p)
ax.set_xlabel("$x (m)$")
ax.set_ylabel("$Presion (kN/m^2)$")

ax = figura.add_subplot(3, 1, 3)
ax.plot(x, d)
ax.set_xlabel("$x (m)$")
ax.set_ylabel("$Densidad (kg/m^3)$")

plt.savefig("estado_"+tiempo+".pdf")
