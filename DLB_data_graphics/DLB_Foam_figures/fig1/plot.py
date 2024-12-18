import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 500,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'lines.markersize': 3,
    'lines.linewidth': 0.7,
    'font.size': 7, # was 10
    'legend.fontsize': 7, # was 10
    'xtick.labelsize': 7,
    'text.usetex': True,
    'ytick.labelsize': 7,
    'figure.figsize': [2.64, 2.5],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)


def generate_target_grid(x,z,xspacing,yspacing):
    xi = np.arange(np.min(x),np.max(x),xspacing)
    zi = np.arange(np.min(z),np.max(z),yspacing)
    xi,zi = np.meshgrid(xi,zi)
    return xi,zi
def interpolate_grid(x,z,field,xi,zi):
    data = griddata((x,z),field,(xi,zi),method='linear')
    return data



cases = ['../../data/shearlayer/scotch/standard_32/']
colors=['b','r']
tmp = open('../../data/shearlayer/scotch/standard_32/processor0/loadBal/cpu_solve.out')
mysize =  np.size(tmp.readlines()[1:])

data = []
for idxcase,case in enumerate(cases):
    mean=0
    rankCounter = 0
    nprocs = np.size(glob.glob(case+'processor*'))
    tmp = []
    for i,rank in enumerate(sorted(glob.glob(case+'processor*'))):
        solve_buffer = []
        time = []
        cpu_solve = open(rank+'/loadBal/cpu_solve.out')
        lines_solve = cpu_solve.readlines()[1:]
        for x in lines_solve:
            time.append(x.split()[0])
            solve_buffer.append(x.split()[4])
        if(i==0):
            size = np.size(solve_buffer)
        solve_buffer = np.array([float(i) for i in solve_buffer])
        tmp.append(solve_buffer[:size-5])
    data.append(tmp)

tmp = np.vstack(data[0])

# Plot 450th datapoint
i = 450
fig,ax = plt.subplots(1)
barlist = ax.bar(np.arange(np.size(tmp[:,0])),tmp[:,i],color='blue')
barlist[np.where(tmp[:,i] == np.max(tmp[:,i]))[0][0]].set_color('red')
ax.axhline(y=np.mean(tmp[:,i]),linestyle='--',color='gray',label='Mean CPU time')
ax.set_xlabel('Processor ID')
ax.set_ylabel('Chemistry CPU time [s]')
ax.set_xlim([0,nprocs])
ax.set_ylim([0,40])
ax.legend(loc=2,frameon=False)
fig.tight_layout()
fig.savefig('chemical_imbalance.pdf',bbox_inches='tight',pad_inches=0.01)
plt.show(block=False)
plt.pause(2)
plt.close("all")


