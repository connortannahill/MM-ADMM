import matplotlib
from string import Template
import subprocess
import glob
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.pyplot import cm
import time
from pathlib import Path
import json

"""
TEMPLATES AND MISC
"""

FUNS = """
create_input()
grid_scale_test_2d()
output_grid_scale_test_2d()
run_parallel_experiment()
create_parallel_plot()
"""

FUNS_LIST = FUNS.split()

def getTemplate2D():
    TEMPLATE_2D = """$testType
    $monType
    $boundaryType
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
    $boundaryType
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

def plot_parallel_experiment(testName, pows, times, ax, color, label):
    
    import statsmodels.stats.api as sms

    # Build each of the confidence intervals
    mean = []
    low_int = []
    high_int = []
    for pow, time_list in times.items():
        interval = sms.DescrStatsW(time_list).tconfint_mean()
        mean.append(np.mean(time_list))
        low_int.append(interval[0])
        high_int.append(interval[1])
    
    ax.plot(np.log2(pows), mean, label=label)
    ax.fill_between(np.log2(pows), low_int, high_int, color=color, alpha=.1)

"""
Functions to be called
"""

def create_parallel_plot():
    testName = input('test name = ')

    # Get all input file names
    dataFilesJson = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/Data/{0}/*'.format(testName))]
    dataFiles = [file[:file.rfind('.')] for file in dataFilesJson]
    num_list = [int(file[len(testName):]) for file in dataFiles]
    sort_list = np.argsort(num_list)
    num_list = list(np.array(num_list)[sort_list])
    dataFiles = list(np.array(dataFiles)[sort_list])
    dataFilesJson = list(np.array(dataFilesJson)[sort_list])

    color=cm.rainbow(np.linspace(0,1,len(dataFiles)))
    num_simplices = [4*(i**2) for i in num_list]
    fig, ax = plt.subplots()

    for i, dataFileJson in enumerate(dataFilesJson):
        times = {}
        print(dataFileJson)
        with open('./Experiments/Data/{0}/{1}'.format(testName, dataFileJson)) as f:
            times = json.loads(f.read())
        
        pows = [int(i) for i in times.keys()]
        
        # Dump the data file
        label = str(num_simplices[i])
        plot_parallel_experiment(dataFiles[i], pows, times, ax, color[i], label)
    
    # Make modified ticks
    ticks = ['${}$'.format(int(2**np.log2(pow))) for pow in pows]
    ax.set_xticklabels(ticks)

    plt.xlabel('Log Number of CPU Cores')
    plt.ylabel('Average CPU Time')
    plt.title('{}'.format(testName))
    plt.legend()

    pName = "Experiments/Results/{0}/".format(testName)
    Path(pName).mkdir(parents=True, exist_ok=True)
    plt.savefig('{0}ParTest{1}.png'.format(pName, testName))
    

def run_parallel_experiment():

    testName = input('test name = ')

    # Get all input file names
    inputFiles = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/InputFiles/{0}*'.format(testName))]
    num_list = np.argsort([int(file[len(testName):]) for file in inputFiles])
    inputFiles = list(np.array(inputFiles)[num_list])[:3]

    HIGHEST_POW = 5

    pows = [2**i for i in range(HIGHEST_POW+1)]

    subprocess.run('make')

    color=cm.rainbow(np.linspace(0,1,len(inputFiles)))
    num_simplices = [4*((2**i * 10)**2) for i in range(1, HIGHEST_POW+5)]
    fig, ax = plt.subplots()
    for i, inputFile in enumerate(inputFiles):
        times = {}
        num_runs = 10
        for pow in pows:
            times[pow] = []
            for run in range(num_runs):
                start = time.time()
                subprocess.run('./mesh.exe {0} {1}'.format(inputFile, pow).split())
                times[pow].append(time.time() - start)
        
        # Dump the data file
        Path("Experiments/Data/{0}/".format(testName)).mkdir(parents=True, exist_ok=True)

        with open('Experiments/Data/{0}/{1}.json'.format(testName, inputFile), 'w+') as f:
            f.write(json.dumps(times))
        
        label = str(num_simplices[i])
        plot_parallel_experiment(inputFile, pows, times, ax, color[i], label)

    plt.xlabel('Number of CPU Cores')
    plt.ylabel('Average Execution times')
    plt.title('{}'.format(testName))
    plt.legend()

    plt.savefig('ParTest{0}.png'.format(testName))
    
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