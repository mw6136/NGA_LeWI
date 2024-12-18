import matplotlib.pyplot as plt
import numpy as np
import matplotlib
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
    'text.usetex': True,
    'ytick.labelsize': 7,
    'figure.figsize': [2.64, 2.4],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)




decomposition = ['simple_x/','simple_y/','scotch/']
label = ['Simple-X','Simple-Y','Scotch']
width = 0.2  # the width of the bars

style = ['-','--','dotted']
cases_loadbal = ['loadbal_32/log.reactingFoam']
cases_noref = ['loadbal_noref_32/log.reactingFoam']
cases_standard = ['standard_32/log.reactingFoam']
x = np.arange(len(label))  # the label locations

fig,ax = plt.subplots(1)
Y1 = []
tmp1 = []
tmp2 = []
tmp3 = []

for decompose in decomposition:
    i = 0
    
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
        tmp1.append(Sr1[-1])
        tmp2.append(Sr2[-1])
        tmp3.append(Sr3[-1])

        counter+=1
    i+=1
rects1 = ax.bar(x-(1.25*width)  , tmp1,width,color='black',fill=False,alpha=0.8, hatch='////',label='LoadBal + Ref',linewidth=0.8)
rects3 = ax.bar(x , tmp3, width,color='black',alpha=0.8,fill=False, label='LoadBal',linewidth=0.8)
rects2 = ax.bar(x +(1.25*width), tmp2, width,color='black',alpha=0.8,hatch='....',fill=False, label='Standard',linewidth=0.8)

ax.legend(loc=2,frameon=False)
ax.set_xticks([0,1,2])
ax.set_xticklabels(['Simple-X','Simple-Y','Scotch'])
ax.set_ylabel('Execution Time [s]')
ax.ticklabel_format(axis='y',style='sci',scilimits=(2,4))
ax.set_ylim([0,None])
ax.set_xlabel('Decomposition Method')
fig.tight_layout()

fig.savefig('balancing_bar.pdf',bbox_inches='tight',pad_inches=0)
plt.show(block=False)
plt.pause(2)
plt.close()