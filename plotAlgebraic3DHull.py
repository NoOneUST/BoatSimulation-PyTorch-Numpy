from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


x, y, L, B, D = 1, 1, 1, 1, 4

Precision = 101
MidIndex = (Precision - 1) / 2

xArray = np.linspace(-0.5, 0.5, Precision)
yArray = np.linspace(-1, 1, Precision)

xGrid, yGrid = np.meshgrid(xArray, yArray)
zGrid = D * (np.power(2 * xGrid / L, 4) + np.power(2 * yGrid / B, 2))
DPlane = np.zeros(zGrid.shape) + D

fig1 = plt.figure(1)
ax = Axes3D(fig1)

ax.plot_surface(xGrid, yGrid, zGrid, rstride=1, cstride=1, cmap='rainbow')
#ax.plot_surface(xGrid, yGrid, DPlane, rstride=1, cstride=1, cmap='rainbow')

plt.show()

fig2 = plt.figure(2)

xNow = xGrid[:, MidIndex]
yNow = yGrid[:, MidIndex]
zNow = zGrid[:, MidIndex]

plt.plot(yNow, zNow)

# plt.show()

#fig3 = plt.figure(3)

xNow = xGrid[:, 0]
yNow = yGrid[:, 0]
zNow = zGrid[:, 0]

plt.plot(yNow, zNow)

plt.show()

fig4 = plt.figure(4)

xNow = xGrid[MidIndex, :]
yNow = yGrid[MidIndex, :]
zNow = zGrid[MidIndex, :]

plt.plot(xNow, zNow)

plt.show()

fig5 = plt.figure(5)

contour = plt.contour(xGrid, yGrid, zGrid, [D], colors='k')

plt.clabel(contour, fontsize=10, colors='k')

#plt.contourf(xGrid, yGrid, zGrid)

#plt.contour(xGrid, yGrid, zGrid)

plt.show()
