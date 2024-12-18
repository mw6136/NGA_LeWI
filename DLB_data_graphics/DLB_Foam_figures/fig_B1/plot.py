import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib

params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 1000,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'font.size': 7, # was 10
    'legend.fontsize': 7, # was 10
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'figure.figsize': [3, 4],
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)

fig,ax = plt.subplots(2,sharex=True)
mean = 0
size = 100
for rankid,rank in enumerate(sorted(glob.glob('../../data/3d_spray/balance/processor*'))):
    get_problem = []
    update_state = []
    balance    = []
    solve_buffer = []
    unbalance = []

    time = []
    cpu_solve = open(rank+'/loadBal/cpu_solve.out')
    lines_solve = cpu_solve.readlines()[1:]
    for x in lines_solve:
        time.append(x.split()[0])
        get_problem.append(x.split()[1])
        update_state.append(x.split()[2])
        balance.append(x.split()[3])
        solve_buffer.append(x.split()[4])
        unbalance.append(x.split()[5])
    time = np.array([float(i) for i in time])
    get_problem = np.array([float(i) for i in get_problem])
    update_state = np.array([float(i) for i in update_state])
    balance = np.array([float(i) for i in balance])
    solve_buffer = np.array([float(i) for i in solve_buffer])
    unbalance = np.array([float(i) for i in unbalance])         
    total = get_problem[1:size]  + update_state[1:size] + balance[1:size] + solve_buffer[1:size] + unbalance[1:size]
    mean += solve_buffer[1:size]
    ax[0].plot(solve_buffer[1:size],linewidth=0.6,label='Processor ID' if rankid==0 else "")
    ax[1].plot((get_problem[1:size]+update_state[1:size]+balance[1:size])*100/(get_problem[1:size]+update_state[1:size]+balance[1:size]+ solve_buffer[1:size]+unbalance[1:size]),linewidth=0.6,label='Processor ID' if rankid==0 else "")

ax[0].plot(mean/np.size(glob.glob('../../data/3d_spray/balance/processor*')),'k--',label='Mean',linewidth=0.9)
ax[0].tick_params(bottom=True,top=False,left=True,right=True,labeltop=False,labelright=False,length=2,direction='in')
ax[0].set_ylabel('Chemistry CPU time [s]')
ax[1].set_ylabel('Overhead [$\%$]')
ax[0].set_ylim([0,25])
ax[0].set_xlim([0,100])
ax[0].legend(loc=1,frameon=False)
fig.align_ylabels(ax[:])
fig.text(0.5, -0.01, 'Number of iterations', ha='center')
fig.tight_layout()
fig.savefig('overhead.pdf',bbox_inches='tight',pad_inches=0)
plt.show(block=False)
plt.pause(2)
plt.close()

