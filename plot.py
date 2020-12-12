import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

points = np.genfromtxt('points.txt', delimiter=',')
triangles = np.genfromtxt('triangles.txt', delimiter=',')
# X = np.genfromtxt('pointsPerfect.txt', delimiter=',')

triang = mtri.Triangulation(points[:,0], points[:,1], triangles=triangles)
# plt.triplot(triang, color='r', marker='o')
plt.triplot(triang, color='r', linewidth=0.5)
# X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
# plt.quiver(X[:,0], X[:,1], points[:,0], points[:,1])

plt.show()
