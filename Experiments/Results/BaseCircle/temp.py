import numpy as np
import glob
files = glob.glob('CircleEx*triangles.txt')
def modCSV(fname):
    vals = np.genfromtxt(fname, delimiter=',', dtype=int)
    vals -= 1
    return vals

for file in files:
    vals = modCSV(file)
    np.savetxt(file, vals, delimiter=',', fmt='%i')