import matplotlib
from string import Template
import subprocess
import glob
import ast
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.pyplot import cm
import time
from pathlib import Path
import json
from os import path

"""
TEMPLATES AND MISC
"""

FUNS = """
create_input()
grid_scale_test_2d()
grid_scale_test_3d()
output_grid_scale_test_2d()
run_parallel_experiment()
run_scale_experiment()
run_simultaneous_experiment()
create_parallel_plot()
plot_single_thread_increase()
plot_energy_decrease()
exit()
"""

FUNS_LIST = FUNS.split()

def getTemplate2D():
    TEMPLATE_2D = """{"TestType": $testType,
    "Dim": $dim,
    "MonType": $monType,
    "CompMesh": $compMesh,
    "BoundaryType": $boundaryType,
    "nSteps": $nSteps,
    "dt": $dt,
    "tau": $tau,
    "rho": $rho,
    "w": $w,
    "nx": $nx,
    "ny": $ny,
    "xa": $xa,
    "xb": $xb,
    "ya": $ya,
    "yb": $yb}"""
    return TEMPLATE_2D

def getTemplate3D():
    TEMPLATE_3D = """{"TestType": $testType,
    "Dim": $dim,
    "MonType": $monType,
    "CompMesh": $compMesh,
    "BoundaryType": $boundaryType,
    "nSteps": $nSteps,
    "dt": $dt,
    "tau": $tau,
    "rho": $rho,
    "w": $w,
    "nx": $nx,
    "ny": $ny,
    "nz": $nz,
    "xa": $xa,
    "xb": $xb,
    "ya": $ya,
    "yb": $yb}
    "za": $za,
    "zb": $zb}"""
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

    plt.savefig("{0}/{1}.png".format(outDir, testName))
    plt.clf()

def exit():
    import sys
    sys.exit()

def plot_parallel_experiment(dirName, testName, num, pows, times, ax, color, mark, label, singlePlot):
    
    import statsmodels.stats.api as sms

    # Build each of the confidence intervals
    mean = []
    low_int = []
    high_int = []
    mean_max = 0
    low_max = 0
    high_max = 0
    i = 0

    tps_mean = []
    tps_low_int = []
    tps_high_int = []
    tps_mean_max = 0
    tps_low_max = 0
    tps_high_max = 0

    denom = 1
    if singlePlot:
        denom = np.mean(times['1'])
    
    # print(times)
    # assert(False)

    for pow, time_list in times.items():
        time_list = [time/denom for time in time_list]
        interval = sms.DescrStatsW(time_list).tconfint_mean()
        mean.append(np.mean(time_list))
        # low_int.append(interval[0])
        # high_int.append(interval[1])
        if i == 0:
            mean_max = mean[0]
            # low_max = low_int[0]
            # high_max = high_int[0]
        i += 1

        # Get num time steps
        nSteps = float(np.genfromtxt('./Experiments/Results/{0}/IhPara{1}.txt'.format(dirName, pow), delimiter=',').shape[0])

        time_mean_list = [time / nSteps for time in time_list]
        # interval = sms.DescrStatsW(time_mean_list).tconfint_mean()
        tps_mean.append(np.mean(time_mean_list))
        # tps_low_int.append(interval[0])
        # tps_high_int.append(interval[1])
        if i == 0:
            tps_mean_max = tps_mean[0]
            # tps_low_max = tps_low_int[0]
            # tps_high_max = tps_high_int[0]
            
    
    if not singlePlot:
        Ih1 = np.genfromtxt('./Experiments/Results/{0}/Ih1.txt'.format(dirName), delimiter=',')
        print('Reading from {}'.format('./Experiments/Results/{0}/Ih1.txt'.format(dirName)))
        eulerTime = (Ih1[-1,:][0])
        print(Ih1)
        # assert(False)

        Ih2 = np.genfromtxt('./Experiments/Results/{0}/Ih2.txt'.format(dirName), delimiter=',')
        BETime = (Ih2[-1,:][0])
        # print(eulerTime)
        # print(BETime)

        ax[0].plot(np.log2(pows), [eulerTime for i in range(len(pows))], label='Euler\'s Method', color='tab:orange')
        ax[0].plot(np.log2(pows), [BETime for i in range(len(pows))], label='Backwards Euler\'s Method', color='g')

        ax[0].plot(np.log2(pows), mean, label=label, color=color, marker=mark)
        # ax[0].fill_between(np.log2(pows), low_int, high_int, color=color, alpha=.1)

        ax[0].legend(prop={'size': 8})

        # ax[1].plot(np.log2(pows), tps_mean, label=label, color=color, marker=mark)
        ax[1].plot(np.log2(pows), tps_mean, color=color, marker=mark)
        # ax[1].fill_between(np.log2(pows), tps_low_int, tps_high_int, color=color, alpha=.1)
        
        # ax[1].legend(prop={'size': 8})
    else:
        print('hit')
        ax[0].plot(np.log2(pows), mean / mean_max, label=label, color=color, marker=mark)#, ms=mark)
        # ax[0].fill_between(np.log2(pows), low_int / low_max, high_int / high_max, color=color, alpha=.1)
        ax[0].legend(prop={'size': 8})

        tps_max = max(tps_mean)
        # ax[1].plot(np.log2(pows), [tps/tps_max for tps in tps_mean], label=label, color=color, marker=mark)
        ax[1].plot(np.log2(pows), [tps/tps_max for tps in tps_mean], color=color, marker=mark)
        # ax[1].legend(prop={'size': 8})

        # ax.plot(np.log2(pows), [tps_mean[0] for i in range(len(pows))] , label=label, color=color)
        # ax[1].fill_between(np.log2(pows), tps_low_int / tps_low_max, tps_high_int / tps_high_max, color=color, alpha=.1)
        # ax.legend(prop={'size': 8})
    # ax.plot(np.log2(pows), mean, label=label, color=color)
    # ax.fill_between(np.log2(pows), low_int, high_int, color=color, alpha=.1)

    return denom

"""
Functions to be called
"""

def plot_energy_decrease():
    testName = input('test name = ')
    timePlot = input('time plot? (True False) ')
    allPlot = input('Plot all experiments? (True False) ')
    print(allPlot)
    methodNum = 0
    if allPlot == 'False':
        methodNum = int(input('Which method? (0 1 2) '))


    vals = []
    if allPlot == 'True':
        for methodType in range(3):
            methodName = ''
            if methodType == 0:
                methodName = 'ADMM'
            elif methodType == 1:
                methodName = 'Euler'
            else:
                methodName = 'Backwards Euler'

            print('./Experiments/Results/{0}/Ih{1}.txt'.format(testName, methodType))
            out =  np.genfromtxt('./Experiments/Results/{0}/Ih{1}.txt'.format(testName, methodType), delimiter=',')
            assert(out.size > 0)
            tVals = out[:,0][1:]
            Ih = out[:,1][1:]
            # with open('./Experiments/Results/{0}/Ih{1}.txt'.format(testName, methodType)) as f:
                # vals = np.array([float(i) for i in f.read().split()])
            if timePlot:
                plt.plot(tVals, Ih, marker='o', ms=2, label='Method Type {}'.format(methodName))
            else:
                plt.plot(np.arange(len(vals))+1, Ih, marker='o', ms=2, label='Method Type {}'.format(methodType))
    else:
        out =  np.genfromtxt('./Experiments/Results/{0}/Ih{1}.txt'.format(testName, methodNum), delimiter=',')
        tVals = out[:,0]
        Ih = out[:,1]

        methodName = ''
        if methodNum == 0:
            methodName = 'ADMM'
        elif methodNum == 1:
            methodName = 'Euler'
        else:
            methodName = 'Backwards Euler'

        # with open('./Experiments/Results/{0}/Ih{1}.txt'.format(testName, methodType)) as f:
            # vals = np.array([float(i) for i in f.read().split()])
        if (timePlot):
            plt.plot(tVals, Ih, marker='o', ms=0.5, label='Method Type {}'.format(methodName))
        else:
            plt.plot(np.arange(len(vals))+1, Ih, marker='o', ms=0.5, label='Method Type {}'.format(methodNum))

    if timePlot == 'True':
        plt.xlabel('CPU time')
    else:
        plt.xlabel('number of time steps')

    NS =  np.genfromtxt('./Experiments/Results/{0}/triangles.txt'.format(testName), delimiter=',').shape[0]
    plt.ylabel('$I_h$')

    mon_name = 'ERR'

    print(testName)

    if testName.startswith('Monitor1'):
        mon_name = 'Example 7.1 with {0} simplices'.format(NS)
    elif testName.startswith('Monitor2'):
        mon_name = 'Example 7.2 with {0} simplices'.format(NS)
    elif testName.startswith('Monitor3'):
        mon_name = 'Example 7.3 with {0} simplices'.format(NS)

    plt.title('{}'.format(mon_name))
    plt.legend()
    plt.savefig('./Experiments/Results/{0}/IhPlot{0}.png'.format(testName))
    # plt.show()

def create_parallel_plot():
    testName = input('test name = ')

    plotAll = input('All on one plot? ') == 'True'
    print(plotAll)

    # Get all input file names
    print('path = {}'.format('./Experiments/Data/{0}/{1}*'.format(testName, testName)))
    print(glob.glob('./Experiments/Data/{0}/{1}*'.format(testName, testName)))
    dataFilesPre = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/Data/{0}/{1}*'.format(testName, testName))]
    dataFilesJson = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/Data/{0}/Para{1}*'.format(testName, testName))]
    dataFiles = [file[:file.rfind('.')] for file in dataFilesJson]
    dataFilesPre = [file[:file.rfind('.')] for file in dataFilesPre]
    num_list = [int(file[len('Para'+testName):]) for file in dataFiles]
    sort_list = np.argsort(num_list)
    num_list = list(np.array(num_list)[sort_list])
    dataFiles = list(np.array(dataFiles)[sort_list])

    dataFilesJson = list(np.array(dataFilesJson)[sort_list])

    num_list2 = [int(file[len(testName):]) for file in dataFilesPre]
    sort_list2 = np.argsort(num_list2)
    dataFilesPre = list(np.array(dataFilesPre)[sort_list2])


    # print(dataFilesPre)
    # print(dataFilesJson)
    # assert(False)


    color=cm.rainbow(np.linspace(0,1,len(dataFiles)))
    marks = ['.', 'o', 'v', '^', 's', 'P']
    # if dim == 2:
    #     num_simplices = [4*(i**2) for i in num_list]
    # else:
    #     num_simplices = [12*(i**3) for i in num_list]
    # if plotAll:
    fig, ax = plt.subplots(2)

    pows = [1, 2, 4, 8, 16]

    for i, dataFileJson in enumerate(dataFilesJson):
        # if not plotAll:
        #     fig, ax = plt.subplots(1 if plotAll else 2)

        times = {}
        print(dataFileJson)
        # with open('./Experiments/Data/{0}/{1}'.format(testName, dataFileJson)) as f:
        #     times = json.loads(f.read())
        
        # pows = [int(i) for i in times.keys()]
        times = {}
        for pow in pows:
            fname = './Experiments/Results/{0}{1}/IhPara{2}.txt'.format(testName, num_list[i], pow)
            if (path.exists(fname)):
                val = np.genfromtxt(fname, delimiter=',')[-1, 0]
                times[str(pow)] = [val]
        print('pows = {}'.format(pows))
        print(times)
        
        # Dump the data file
        # label = str(num_simplices[i])
        # Get the number of simplces
        num_simplices = np.genfromtxt('./Experiments/Results/{0}{1}/triangles.txt'.format(testName, num_list[i])).shape[0]
        label = 'MM-ADMM on ' + str(num_simplices) + ' simplices'
        print('label {}'.format(label))
        denom = plot_parallel_experiment(dataFilesPre[i], dataFiles[i], num_list[i], pows, times, ax, color[i], marks[i], label, plotAll)

        if not plotAll:
            # Make modified ticks
            # ticks = ['${}$'.format(pow) for pow in pows]
            # ticks.insert(0, '$1$')
            # ax[0].set_xticklabels(ticks)
            # ax[1].set_xticklabels(ticks)
            plt.setp(ax, xticks=[i for i in range(len(pows))], xticklabels=[str(pow) for pow in pows])

            plt.xlabel('Number of CPU Cores')
            # plt.ylabel('Average CPU Time')
            ax[0].set(ylabel='Average CPU Time')
            ax[1].set(ylabel='Average CPU TPS')
            # ax[0].ylabel('Average CPU Time')
            # ax[1].ylabel('Average CPU Time Per Step')
            # plt.title('{}'.format(testName))
            # plt.legend()

            pName = "Experiments/Results/{0}{1}/".format(testName, num_list[i])
            print(pName)
            Path(pName).mkdir(parents=True, exist_ok=True)
            plt.savefig('{0}ParTest{1}{2}.png'.format(pName, testName, num_list[i]))
            plt.clf()

    if plotAll:
        # ticks = ['${}$'.format(pow) for pow in pows]
        # ticks.insert(0, '$1$')
        # print(ticks)
        # ax.set_xticklabels(ticks)
        # plt.xlabel('Number of CPU Cores')
        # plt.ylabel('Normalized Average CPU Time')
        # plt.title('{}'.format(testName))
        # plt.legend()

        # pName = "Experiments/Results/"
        # Path(pName).mkdir(parents=True, exist_ok=True)
        # plt.savefig('{0}ParTest{1}.png'.format(pName, testName))

        # ticks = ['${}$'.format(pow) for pow in pows]
        # ticks.insert(0, '$1$')
        # print(ticks)
        # ax[0].ticks(np.arange(len(pow)))
        # ax[0].set_xticklabels(ticks)
        # ax[1].ticks(np.arange(len(pow)))
        # ax[1].set_xticklabels(ticks)

        plt.setp(ax, xticks=[i for i in range(len(pows))], xticklabels=[str(pow) for pow in pows])

        plt.xlabel('Number of CPU Cores')
        ax[0].set(ylabel='Normalized CPU Time')
        ax[1].set(ylabel='Normalized CPU TPS')
        # plt.ylabel('Normalized Average CPU Time')
        # plt.title('{}'.format(testName))
        if testName.startswith('Monitor1'):
            mon_name = 'Example 7.1 Parallel Scaling'
        elif testName.startswith('Monitor2'):
            mon_name = 'Example 7.2 Parallel Scaling'
        elif testName.startswith('Monitor3'):
            mon_name = 'Example 7.3 Parallel Scaling'

        fig.suptitle('{}'.format(mon_name))
        # fig.title()
        # plt.legend()

        pName = "Experiments/Results/"
        Path(pName).mkdir(parents=True, exist_ok=True)
        plt.savefig('{0}ParTest{1}.png'.format(pName, testName))
    
def run_parallel_experiment():
    HIGHEST_POW = 5

    testName = input('test name = ')

    # Get all input file names
    inputFiles = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/InputFiles/{0}*'.format(testName))]
    num_literal = [int(file[len(testName):file.rfind('.')]) for file in inputFiles]
    num_list = np.argsort(num_literal)
    num_literal = np.sort(num_literal)
    inputFiles = [s[:s.rfind('.')] for s in list(np.array(inputFiles)[num_list])]

    subprocess.run('make')

    pows = [2**i for i in range(HIGHEST_POW+1)]

    for i, inputFile in enumerate(inputFiles):
        times = {}
        num_runs = 10

        for pow in pows:
            times[pow] = []
            for run in range(num_runs):
                start = time.time()
                subprocess.run('./mesh.exe {0} 0 {1}'.format(inputFile, pow).split())
                Ih = np.genfromtxt('./Experiments/Results/{0}/Ih0.txt'.format(inputFile), delimiter=',')
                runTime = (Ih[-1,:][0])
                times[pow].append(runTime)
        
        # Dump the data file
        Path("Experiments/Data/{0}/".format(testName)).mkdir(parents=True, exist_ok=True)

        with open('Experiments/Data/{0}/Para{1}.json'.format(testName, inputFile), 'w+') as f:
            f.write(json.dumps(times))

def run_simultaneous_experiment():

    HIGHEST_POW = 5

    testName = input('test name = ')

    # Get all input file names
    inputFiles = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/InputFiles/{0}*'.format(testName))]
    num_literal = [int(file[len(testName):file.rfind('.')]) for file in inputFiles]
    num_list = np.argsort(num_literal)
    num_literal = np.sort(num_literal)
    inputFiles = [s[:s.rfind('.')] for s in list(np.array(inputFiles)[num_list])]

    subprocess.run('make')

    pows = [2**i for i in range(HIGHEST_POW+1)]

    for i, inputFile in enumerate(inputFiles):
        times = {}
        num_runs = 10

        times['({0}, {1})'.format(num_list[i], pows[i])] = []
        for run in range(num_runs):
            start = time.time()
            subprocess.run('./mesh.exe {0} 0 {1}'.format(inputFile, pows[i]).split())
            times['({0}, {1})'.format(num_list[i], pows[i])].append(time.time() - start)
        
        # Dump the data file
        Path("Experiments/Data/{0}/".format(testName)).mkdir(parents=True, exist_ok=True)

        with open('Experiments/Data/{0}/Simul{1}.json'.format(testName, inputFile), 'w+') as f:
            f.write(json.dumps(times))

def run_scale_experiment():

    testName = input('test name = ')

    # Get all input file names
    inputFiles = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/InputFiles/{0}*'.format(testName))]
    print(inputFiles)
    num_literal = [int(file[len(testName):file.rfind('.')]) for file in inputFiles]
    print(num_literal)
    num_list = np.argsort(num_literal)
    print(num_list)
    num_literal = np.sort(num_literal)
    print(num_list)
    inputFiles = [s[:s.rfind('.')] for s in list(np.array(inputFiles)[num_list])]

    # HIGHEST_POW = 5

    # pows = [2**i for i in range(HIGHEST_POW+1)]

    subprocess.run('make')

    fig, ax = plt.subplots()
    print(inputFiles)
    for i, inputFile in enumerate(inputFiles):
        times = {i: {} for i in range(3)}
        for method in range(3):
            num_runs = 1
            # for pow in pows:
            times[method][int(num_literal[i])] = []
            for run in range(num_runs):
                start = time.time()
                subprocess.run('./mesh.exe {0} {1}'.format(inputFile, method).split())
                times[method][int(num_literal[i])].append(time.time() - start)
            
        # Dump the data file
        Path("Experiments/Data/{0}/".format(testName)).mkdir(parents=True, exist_ok=True)

        with open('Experiments/Data/{0}/Single{1}.json'.format(testName, inputFile), 'w+') as f:
            f.write(json.dumps(times))
        
def plot_single_thread_increase():

    testName = input('test name = ')
    
    color = cm.rainbow(np.linspace(0,1,5))


    # Add the single thread average time for each experiment to the plot
    # Get all of the data files

    print('path = {}'.format('./Experiments/Data/{0}/{1}*'.format(testName, testName)))
    print(glob.glob('./Experiments/Data/{0}/{1}*'.format(testName, testName)))
    dataFilesJson = [file[file.rfind('/')+1:] for file in glob.glob('./Experiments/Data/{0}/{1}*'.format(testName, testName))]
    dataFiles = [file[:file.rfind('.')] for file in dataFilesJson]
    num_list = [int(file[len(testName):]) for file in dataFiles]
    sort_list = np.argsort(num_list)
    num_list = list(np.array(num_list)[sort_list])
    dataFiles = list(np.array(dataFiles)[sort_list])
    dataFilesJson = list(np.array(dataFilesJson)[sort_list])

    times_list_0 = []
    times_list_1 = []
    times_list_2 = []
    simplices = []
    for i, dataFileJson in enumerate(dataFilesJson):
        num_simplices = (np.genfromtxt('./Experiments/Results/{0}{1}/triangles.txt'.format(testName, num_list[i]))).shape[0]
        simplices.append(num_simplices)
        times = {}
        with open('./Experiments/Data/{0}/{1}'.format(testName, dataFileJson)) as f:
            times = json.loads(f.read())
        
        times_list_0.append(np.mean(times['0'][str(num_list[i])]))
        times_list_1.append(np.mean(times['1'][str(num_list[i])]))
        times_list_2.append(np.mean(times['2'][str(num_list[i])]))
        
        # one_t_time = times[methodName][str(num_list[i])]
        # print(one_t_time)
        # times_list.append(np.mean(one_t_time))
    

    # Scatter plot of the times
    mask = np.array(num_list) < 320
    num_list = np.array(num_list)[mask]
    simplices = np.array(simplices)[mask]
    times_list_0 = np.array(times_list_0)[mask]
    times_list_1 = np.array(times_list_1)[mask]
    times_list_2 = np.array(times_list_2)[mask]

    # times_list = np.array(times_list)[num_list<63]
    plt.plot(simplices, times_list_0, color=color[0], label='MM-ADMM')
    plt.plot(simplices, times_list_1, color=color[1], label='Euler\'s Method')
    plt.plot(simplices, times_list_2, color=color[2], label='Backward Euler\'s Method')

    # # Build regression
    # A = np.vstack([num_list**3, num_list**2, num_list, np.ones(len(num_list))]).T
    # coefs, res, rank, s = np.linalg.lstsq(A, times_list, rcond=None)
    # # a, b, c = coefs
    # a, b, c, d = coefs
    # x = np.linspace(min(num_list), max(num_list), 500)
    # plt.plot(x, a*(x**3) + b*(x**2) + c*x + d, color=color[j], label='{0} ($SSE = {1:.2e}, a = {2:.2e}, b = {3:.2e}$)'.format(testName, float(res), float(a), float(b)))

    plt.legend()

    plt.savefig('./Experiments/Results/{0}SingleThread.png'.format(testName))

def output_grid_scale_test_2d():
    testName = input('test name = ')

    out_dirs = glob.glob('./Experiments/Results/*{}*'.format(testName))

    for dir in out_dirs:
        outputPng(dir)

def grid_scale_test_2d():
    paramList = [s[1:-1] for s in getTemplate2D().split()[1::2]]
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
            OUT_DIR + outFileName + str(n) + '.json')

        in_files.append(outFileName + str(n))
    
def grid_scale_test_3d():
    paramList = [s[1:-1] for s in getTemplate3D().split()[1::2] if s[0] == '$']
    paramList.remove('nx')
    paramList.remove('ny')
    paramList.remove('nz')

    in_dict = {s: input('{} = '.format(s)) for s in paramList}

    nVals = [str((2**i)*10) for i in range(6)]

    print(nVals)

    outFileName = input('test name = ')

    in_files = []

    for n in nVals:
        in_dict['nx'] = n
        in_dict['ny'] = n
        in_dict['nz'] = n

        create_input_from_dict(in_dict, getTemplate3D(),
            OUT_DIR + outFileName + str(n))

        in_files.append(outFileName + str(n))

def create_input():
    in_dict = {}

    D = input('DIM = ')

    paramList = [s[1:-1] for s in getTemplate2D().split()[1::2]]
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