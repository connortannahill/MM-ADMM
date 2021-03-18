import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys

mode = int(sys.argv[1])

if mode == 0:

    points = np.genfromtxt('points.txt', delimiter=',')
    # print(points)
    triangles = np.genfromtxt('triangles.txt', delimiter=',')
    # X = np.genfromtxt('pointsPerfect.txt', delimiter=',')

    triang = mtri.Triangulation(points[:,0], points[:,1], triangles=triangles)
    # # plt.triplot(triang, color='r', marker='o')
    plt.triplot(triang, color='r', linewidth=0.5)
    # plt.scatter(points[:,0], points[:,0])
    # X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
    # plt.quiver(X[:,0], X[:,1], points[:,0], points[:,1])

    # import plotly.express as px
    # fig = px.scatter(x=points[:,0], y=points[:,1])
    # fig.show()
elif mode == 1:
    centers = np.genfromtxt('centroid.txt', delimiter=',')

    x = centers[:,0]
    y = centers[:,1]

    print(x)

    plt.scatter(x, y, s=1)
elif mode == 2:
    from mpl_toolkits.mplot3d import Axes3D 

    centers = np.genfromtxt('centroid.txt', delimiter=',')

    x = centers[:,0]
    y = centers[:,1]
    z = centers[:,2]

    print(z)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z)





plt.show()
