import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import pdb
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 1000,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'lines.markersize': 3,
    'lines.linewidth': 0.7,
    'font.size': 7, # was 10
    'legend.fontsize': 7, # was 10
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'text.usetex': True,
    'figure.figsize': [2.64, 2.9],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)


decomposition = ['simple_x/','simple_y/','scotch/']
label = ['Simple-X','Simple-Y','Scotch']
style = ['-','--','dotted']

cases_loadbal = ['loadbal_256/log.reactingFoam','loadbal_128/log.reactingFoam','loadbal_64/log.reactingFoam','loadbal_32/log.reactingFoam']
cases_noref = ['loadbal_noref_256/log.reactingFoam','loadbal_noref_128/log.reactingFoam','loadbal_noref_64/log.reactingFoam','loadbal_noref_32/log.reactingFoam']
cases_standard = ['standard_256/log.reactingFoam','standard_128/log.reactingFoam','standard_64/log.reactingFoam','standard_32/log.reactingFoam']




fig,ax = plt.subplots(1)

for i,decompose in zip(np.arange(np.size(decomposition)),decomposition):
    counter = 0
    Sr1 = []
    Sr2 = []
    Sr3 = []

    for file1,file2,file3 in zip(cases_loadbal,cases_standard,cases_noref):
        files_loadbal = open('../../data/shearlayer/'+decompose+file1,'r')
        files_standard = open('../../data/shearlayer/'+decompose+file2,'r')
        files_noref = open('../../data/shearlayer/'+decompose+file3,'r')

        string = 'ClockTime'
        list1 = []
        list2 = [0]
        for line in files_loadbal:
            if string in line:
                list1.append(float(line.split(' ')[2]))
                list2.append(float(line.split(' ')[2]))    
        Sr1.append((list1[-1]))
        list1 = []
        list2 = [0]
        for line in files_standard:
            if string in line:
                list1.append(float(line.split(' ')[2]))
                list2.append(float(line.split(' ')[2]))
        Sr2.append((list1[-1]))

        list1 = []
        list2 = [0]
        for line in files_noref:
            if string in line:
                list1.append(float(line.split(' ')[2]))
                list2.append(float(line.split(' ')[2]))
        Sr3.append((list1[-1]))

        counter+=1
    Sr1 = np.array(Sr1[::-1])
    Sr2 = np.array(Sr2[::-1])
    Sr3 = np.array(Sr3[::-1])
    ax.plot(Sr2/Sr3,'ks',linestyle=style[i],markerfacecolor='None',markeredgewidth=0.7,linewidth=0.7)
    ax.plot(Sr2/Sr1,'ro',linestyle=style[i],markerfacecolor='None',markeredgewidth=0.7,linewidth=0.7)

    ax.plot(np.NaN, np.NaN, color='black', marker=None,linestyle=style[i], label=label[i])

    ax.set_xlim([-0.1,3.1])
    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels(['32','64','128','256'])

    ax.set_xlabel(r'$N_{procs}$')
    ax.set_ylabel(r'$\chi_{su}$')
ax.plot(np.NaN, np.NaN, color='black', marker='s',markerfacecolor='None',markeredgewidth=0.7,linestyle=None, linewidth=0,label='LoadBal')
ax.plot(np.NaN, np.NaN, color='red', marker='o',markerfacecolor='None',markeredgewidth=0.7,linestyle=None, linewidth=0,label='LoadBal + RefMap')

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),
          frameon=False, shadow=False, ncol=2)
fig.tight_layout()
fig.savefig('speedup.pdf',bbox_inches='tight',pad_inches=0)
plt.show(block=False)
plt.pause(2)
plt.close()