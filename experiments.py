import matplotlib as plt
from string import Template
import subprocess
import glob
import matplotlib.pyplot
import matplotlib.tri as mtri
import numpy as np

"""
TEMPLATES AND MISC
"""

FUNS = """
create_input()
grid_scale_test_2d()
output_grid_scale_test_2d()
"""

FUNS_LIST = FUNS.split()

def getTemplate2D():
    TEMPLATE_2D = """$testType
    $monType
    $nSteps
    $dt
    $tau
    $rho
    $nx
    $ny
    $xa
    $xb
    $ya
    $yb"""
    return TEMPLATE_2D

def getTemplate3D():
    TEMPLATE_3D = """$testType
    $monType
    $nSteps
    $dt
    $tau
    $rho
    $nx
    $ny
    $nz
    $xa
    $xb
    $ya
    $yb
    $za
    $zb"""
    return TEMPLATE_3D

OUT_DIR = './Experiments/InputFiles/'

"""
Helper functions
"""
def create_input_from_dict(in_dict, template, outFileName):
    s = Template(template)
    inFile = s.substitute(in_dict)

    with open(outFileName, 'w') as f:
        f.write(inFile)

def run_experiments(exp_list):
    for exp in exp_list:
        subprocess.run("./mesh.exe {}".format(exp).split())

def outputPng(outDir):
    testName = outDir[outDir.rfind('/')+1:]
    points = np.genfromtxt('{}/points.txt'.format(outDir), delimiter=',')
    triangles = np.genfromtxt('{}/triangles.txt'.format(outDir), delimiter=',')

    triang = mtri.Triangulation(points[:,0], points[:,1], triangles=triangles)
    plt.triplot(triang, color='r')

    plt.savefig("./Experiments/PngOutput/{}.png".format(testName))

"""
Functions to be called
"""

def output_grid_scale_test_2d():
    testName = input('test name = ')

    out_dirs = glob.glob('./Experiments/Results/*{}*'.format(testName))

    for dir in out_dirs:
        outputPng(dir)

def grid_scale_test_2d():
    paramList = [s[1:] for s in getTemplate2D().split()]
    paramList.remove('nx')
    paramList.remove('ny')

    in_dict = {s: input('{} = '.format(s)) for s in paramList}

    nVals = [str((2**i)*10) for i in range(6)]

    print(nVals)

    outFileName = input('test name = ')

    in_files = []

    for n in nVals:
        in_dict['nx'] = n
        in_dict['ny'] = n

        create_input_from_dict(in_dict, getTemplate2D(),
            OUT_DIR + outFileName + str(n))

        in_files.append(outFileName + str(n))
    
    run_experiments(in_files)

def create_input():
    in_dict = {}

    D = input('DIM = ')

    paramList = [s[1:] for s in getTemplate2D().split()]
    in_dict = {s: input('{} = '.format(s)) for s in paramList}

    outFileName = input('out file name = ')

    s = Template(getTemplate2D())
    inFile = s.substitute(in_dict)

    with open(OUT_DIR + outFileName, 'w') as f:
        f.write(inFile)

while(True):
    print('====================================')
    print("Functions available:")
    print(FUNS)
    print('====================================')

    mode = input('call function() = ')
    if mode in FUNS_LIST:
        eval(mode)
    else:
        print('Function {} does not exist!'.format(mode))