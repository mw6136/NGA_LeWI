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
    'text.usetex': True,
    'figure.figsize': [2.64, 5],
    'font.family': 'serif',
}
matplotlib.rcParams.update(params)


cases = ['scotch/loadbal_32/','scotch/loadbal_noref_32/','scotch/standard_32/']
cases = ['scotch/standard_32/','scotch/loadbal_noref_32/','scotch/loadbal_32/']

labels = ['Standard','LoadBal','LoadBal + Ref']
fig,ax = plt.subplots(ncols=1,nrows=np.size(cases),sharex=True,sharey=True)

caseCounter = 0

for ind,case in enumerate(cases):
    tmp  = 0
    mean = 0
    rankCounter = 0
    for rankid,rank in enumerate(sorted(glob.glob('../../data/shearlayer/'+case+'processor*'))):    
        get_problem = []
        update_state = []
        balance    = []
        solve_buffer = []
        unbalance = []
        time = []
        cpu_solve = open('../../data/shearlayer/'+rank+'/loadBal/cpu_solve.out')
        lines_solve = cpu_solve.readlines()[1:]
        for x in lines_solve:
            time.append(x.split()[0])
            get_problem.append(x.split()[1])
            update_state.append(x.split()[2])
            balance.append(x.split()[3])
            solve_buffer.append(x.split()[4])
            unbalance.append(x.split()[5])

        if(rankCounter==0):
            size = np.size(time)-2
        rankCounter+=1
        time = np.array([float(i) for i in time])
        get_problem = np.array([float(i) for i in get_problem])
        update_state = np.array([float(i) for i in update_state])
        balance = np.array([float(i) for i in balance])
        solve_buffer = np.array([float(i) for i in solve_buffer])
        unbalance = np.array([float(i) for i in unbalance])         
        total = get_problem[:size]  + update_state[:size] + balance[:size] + solve_buffer[:size] + unbalance[:size]
        mean += solve_buffer[:size]
        ax[ind].plot(solve_buffer[0::2] + solve_buffer[1::2],linewidth=0.6,label='Processor ID' if rankid==0 else "")
    ax[ind].plot((mean[0::2] + mean[1::2])/np.size(glob.glob('../../data/shearlayer/'+case+'processor*')),'k--',label='Mean',linewidth=0.9)
    ax[ind].tick_params(bottom=True,top=False,left=True,right=True,labeltop=False,labelright=False,length=2,direction='in')
    ax[ind].text(0.03,0.96,'{}'.format(labels[ind]), transform = ax[ind].transAxes,fontsize=8,color='black',horizontalalignment='left',verticalalignment='top')
    ax[ind].set_ylabel('Chemistry CPU time [s]')


ax[0].set_xlim([0,None])
ax[0].legend(loc=1,frameon=False)
ax[0].set_xlim([0,1000])
fig.text(0.5, -0.01, 'Number of iterations', ha='center')
fig.tight_layout()
fig.savefig('rankbased_solve.pdf',bbox_inches='tight',pad_inches=0)
plt.show(block=False)
plt.pause(2)
plt.close()
