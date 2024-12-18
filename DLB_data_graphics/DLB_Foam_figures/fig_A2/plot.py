import scipy.stats
import numpy as np
import matplotlib
import pdb
import matplotlib.pyplot as plt
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 1000,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'font.size': 7, # was 10
    'legend.fontsize': 6, # was 10
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'text.usetex': True,
    'figure.figsize': [2.64, 2.2],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)


proc_counts = [64]
cases = ['loadbal','loadbal_noref','standard']

for count in proc_counts:
    tmp_qdot = []
    for case in cases:
        time_qdot = []
        qdot = []
        path = '../../data/appendix/appendixA/scotch_1e4/{}_{}/postProcessing/Integral_Qdot/0/volFieldValue.dat'.format(case,count)
        with open(path) as f:
            next(f)
            next(f)
            next(f)
            next(f)
            for line in f:
                time_qdot.append(float(line.split()[0]))
                qdot.append(float(line.split()[1]))

        tmp_qdot.append(np.array(qdot))

    fig,ax = plt.subplots()
    ax.plot(tmp_qdot[0],color='black',linestyle='-',lw=0.8,label='LoadBal+Ref')
    ax.plot(tmp_qdot[1],color='red',linestyle='--',lw=0.8,label='LoadBal')
    ax.plot(tmp_qdot[2],color='blue',linestyle='dotted',lw=0.8,label='Standard')
    ax.set_xlabel('Number of iterations')
    ax.set_ylabel('$\int_V \dot{Q}_{i}dV$ [J/s]')
    ax.set_xlim([0,500])
    ax.legend(loc='best',frameon=False)
    fig.tight_layout()
    fig.savefig('refmap_err.pdf',bbox_inches='tight',pad_inches=0)
    plt.show(block=False)
    plt.pause(2)
    plt.close()
