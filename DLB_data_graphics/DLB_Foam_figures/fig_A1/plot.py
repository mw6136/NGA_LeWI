from sklearn.utils import check_array
import scipy.stats
import numpy as np
import matplotlib
import random
import pdb
import sys
import random
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'font.size': 7, # was 10
    'legend.fontsize': 6, # was 10
    'xtick.labelsize': 7,
    'text.usetex': True,
    'ytick.labelsize': 7,
    'figure.figsize': [2.64, 2.2],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)

def mean_absolute_percentage_error(y_true, y_pred): 
    difference = np.abs(y_true-y_pred)
    return np.abs(np.mean(difference / (y_true+np.finfo(float).eps))) * 100

def max_absolute_percentage_error(y_true, y_pred):
    difference = np.abs(y_true-y_pred)
    return np.abs(np.max(difference / (y_true+np.finfo(float).eps))) * 100

def plot_error(tolerance,model,standard,ax,color,marker,label):
    proc_counts = [32,64,128,256]

    error_mapped_qdot = []

    for idx,count in enumerate(proc_counts):
        tmp_qdot = []
        time_qdot = []

        qdot = []
        path = '../../data/appendix/appendixA/{}/{}_{}/postProcessing/Integral_Qdot/0/volFieldValue.dat'.format(tolerance,model,count)
        with open(path) as f:
            next(f)
            next(f)
            next(f)
            next(f)
            for line in f:
                time_qdot.append(float(line.split()[0]))
                qdot.append(float(line.split()[1]))

        tmp_qdot.append(np.array(qdot))

        error_mapped_qdot.append(mean_absolute_percentage_error(tmp_qdot[0],standard[idx]))

    ax.plot(error_mapped_qdot,lw=0.8, color=color, marker=marker,markersize=4,label=label)
    ax.legend()

def main():
    fig,ax = plt.subplots()
    proc_counts = [32,64,128,256]
    # Load standard first
    standard = []
    for count in proc_counts:
        tmp_qdot = []
        tmp_T = []
        time_qdot = []

        qdot = []
        path = '../../data/appendix/appendixA/scotch_1e4/standard_{}/postProcessing/Integral_Qdot/0/volFieldValue.dat'.format(count)
        with open(path) as f:
            next(f)
            next(f)
            next(f)
            next(f)
            for line in f:
                time_qdot.append(float(line.split()[0]))
                qdot.append(float(line.split()[1]))

        tmp_qdot.append(np.array(qdot))
        standard.append(tmp_qdot[0])

    plot_error('scotch_1e4','loadbal_noref',standard,ax,'k','x','No mapping')
    plot_error('scotch_1e2','loadbal',standard,ax,'b','o','Z$_{tol} = 1e^{-2}$')
    plot_error('scotch_1e4','loadbal',standard,ax,'r','o','Z$_{tol} = 1e^{-4}$')
    plot_error('scotch_1e6','loadbal',standard,ax,'g','o','Z$_{tol} = 1e^{-6}$')

    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels(['32','64','128','256'])
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.set_xlabel('$N_{p}$')
    ax.set_ylabel(" err($\int_V \dot{Q}_{i}dV$)")
    ax.legend(loc='best',frameon=False)
    fig.tight_layout()
    fig.savefig('error_qdot.pdf',bbox_inches='tight',pad_inches=0.0)
    plt.show(block=False)
    plt.pause(2)
    plt.close()
main()
