
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import sys, copy
mode = int(sys.argv[1])

import plotly.plotly as py
import plotly.graph_objs as go

import matplotlib.cm as cm
from scipy.spatial import Delaunay

if mode == 0:

    testName = sys.argv[2]
    outDir = 'Experiments/Results/' + testName

    # out = np.genfromtxt('{}/phi.txt'.format(outDir), delimiter=',')

    # x = out[:,0]
    # y = out[:,1]
    # phi = out[:,2]

    # n = int(np.sqrt(phi.size))

    # fig, ax = plt.subplots(1)
    # img = ax.contour(np.reshape(x, (n, n)), np.reshape(y, (n, n)), np.reshape(phi, (n, n)), levels=[0], colors='b')

    points = np.genfromtxt('{}/points.txt'.format(outDir), delimiter=',')
    # print(points)
    triangles = np.genfromtxt('{}/triangles.txt'.format(outDir), delimiter=',')
    # X = np.genfromtxt('pointsPerfect.txt', delimiter=',')

    triang = mtri.Triangulation(points[:,0], points[:,1], triangles=triangles)
    # # plt.triplot(triang, color='r', marker='o')
    pntSet = {pnt for pnt in triangles.flatten()}
    pnts = np.array(list(pntSet), dtype=int)
    # print(pntSet)

    plt.triplot(triang, color='r', linewidth=0.1)
    usedPnts = points[pnts,:]
    plt.scatter(usedPnts[:,0], usedPnts[:,1], c='r', s=0.1)
    # X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
    # plt.quiver(X[:,0], X[:,1], points[:,0], points[:,1])



    # import plotly.express as px
    # fig = px.scatter(x=points[:,0], y=points[:,1])
    # fig.show()
    plt.savefig("sample_plot.png")
elif mode == 1:
    gridData = np.genfromtxt('gridMesh.txt', delimiter=',')

    x = gridData[:,0]
    y = gridData[:,1]
    M = gridData[:,2]

    n = int(np.sqrt(x.size))

    M = np.reshape(M, (n, n))

    plt.imshow(M)
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
elif mode == 3:
    points = np.genfromtxt('boundaryPnts.txt', delimiter=',')
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    ax.scatter(points[:,0], points[:,1], points[:,2])

    # import meshio
    # meshio.write('out.stl', points, {'tetra': triangles})
elif mode == 4:
    """ Create a gif of the input data """
    import os, re
    import numpy as np
    import matplotlib.pyplot as plt
    import imageio
    
    outdir = './gifout/'

    files = os.listdir(outdir)

    rx = re.compile('X.*')
    rz = re.compile('Z.*')

    xNames = list(filter(rx.match, files))
    xNames = [outdir + xName for xName in xNames]

    zNames = list(filter(rz.match, files))
    zNames = [outdir + zName for zName in zNames]

    print(xNames)
    print(zNames)

    triangles = np.genfromtxt('triangles.txt', delimiter=',')

    """ Build the image set for the mesh """

    xFiles = []

    # Generating the x plots
    print('genning x')
    for xName in xNames:
        x = np.genfromtxt(xName, delimiter=',')

        triang = mtri.Triangulation(x[:,0], x[:,1], triangles=triangles)

        plt.triplot(triang, color='r', linewidth=0.5)

        nums = [int(i) for i in re.findall(r'\d+', xName)]
        i = nums[0]

        xFiles.append('X{}.png'.format(i))

        plt.savefig(xFiles[-1])

        plt.clf()
    print('done genning x')

    
    # Generating the z plots

    # Building up the connectivity
    length = triangles.shape[0]
    triangles = []
    # triangles[0]
    for i in range(length):
       triangles.append([3*i, 3*i+1, 3*i+2]) 

    triangles = np.asarray(triangles)

    zFiles = []
    for zName in zNames:
        z = np.genfromtxt(zName, delimiter=',')

        triang = mtri.Triangulation(z[:,0], z[:,1], triangles=triangles)

        plt.triplot(triang, color='b', linewidth=0.5)

        nums = [int(i) for i in re.findall(r'\d+', zName)]

        zFiles.append('Z{0}-{1}.png'.format(nums[0], nums[1]))

        plt.savefig(zFiles[-1])

        plt.clf()
    
    # Sort the lists
    atoi = lambda x: int(x) if x.isdigit() else x
    natural_key = lambda x: [atoi(c) for c in re.split(r'(\d+)', x)]
    
    xFiles.sort(key=natural_key)
    zFiles.sort(key=natural_key)

    print(xFiles)
    print(zFiles)

    xFiles_cpy = copy.deepcopy(xFiles)
    zFiles_cpy = copy.deepcopy(zFiles)

    with imageio.get_writer('meshgif.gif', mode='I') as writer:
        # for xFile in xFiles:
        #     image = imageio.imread(xFile)
        #     writer.append_data(image)

        n = len(xFiles)

        xFile = ''
        for i in range(n):
            # Display the x value
            xFile = xFiles.pop(0)
            image = imageio.imread(xFile)
            writer.append_data(image)

            # Extract z sub-list for appropriate step
            
            # while (int(re.findall(r'\d+', zFiles[0])[0]) == i):
            #     # Display the z value
            #     zFile = zFiles.pop(0)
            #     image = imageio.imread(zFile)
            #     writer.append_data(image)

            #     if (len(zFiles) == 0):
            #         break

        
        # Write 30 redundant frames
        for i in range(30):
            image = imageio.imread(xFile)
            writer.append_data(image)

    # Remove files
    for filename in set(xFiles_cpy):
        os.remove(filename)

    for filename in set(zFiles_cpy):
        os.remove(filename)
elif mode == 5:
    testName = sys.argv[2]
    outDir = 'Experiments/Results/' + testName

    points = np.genfromtxt('{}/points.txt'.format(outDir), delimiter=',')
    triangles = np.genfromtxt('{}/triangles.txt'.format(outDir), delimiter=',')

    # x, y, z, i, j, k = [[] for i in range(6)]
    x = list(points[:, 0])
    y = list(points[:, 1])
    z = list(points[:, 2])

    i, j, k = [[] for i in range(3)]


    for tri in triangles:
        id0, id1, id2, id3 = [int(t) for t in tri]

        possibleFaces = [
            [id0, id1, id2],
            [id0, id1, id3],
            [id0, id2, id3],
            [id1, id2, id3]
        ]

        # Each of the face options
        for face in possibleFaces:
            i.append(face[0])
            j.append(face[1])
            k.append(face[2])

    fig = go.Figure(data=[
        go.Mesh3d(
            x=x
            y=x
            z=x
            i = i,
            j = j,
            k = k,
        )
    ])

    fig.show()
if mode != 4:
    plt.show()

# %%
