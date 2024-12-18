import matplotlib.pyplot as plt
import numpy as np
import matplotlib
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 300,  # to adjust notebook inline plot size
    'axes.labelsize': 7, # fontsize for x and y labels (was 10)
    'axes.titlesize': 7,
    'lines.markersize': 3,
    'lines.linewidth': 0.7,
    'font.size': 7, # was 10
    'legend.fontsize': 6, # was 10
    'xtick.labelsize': 7,
    'text.usetex': True,
    'ytick.labelsize': 7,
    'figure.figsize': [2.64, 2.4],
    'font.family': 'serif',
}
matplotlib.rcParams.update(params)


def get_data(file):

    execution_time = []
    ueqn_time = []
    peqn_time = []
    yeqn_time = []
    eeqn_time = []
    lagrangian_time = []

    myfile = open(file,'r')

    # First get total execution time
    execute = 'ClockTime'

    # U Equation time
    ueqn = 'U Equation Time'

    # Y Equation time
    yeqn = 'Y Equation Time'

    # p Equation time
    peqn = 'p Equation Time'

    # E	Equation time
    eeqn = 'E Equation Time'

    # Lagrangian time
    lagrangian = 'Lagrangian Time'

    for line in myfile:
        if execute in line:
            execution_time.append(float(line.split(' ')[2]))
        elif ueqn in line:
            ueqn_time.append(float(line.split()[-1]))
       	elif	yeqn in	line:
            yeqn_time.append(float(line.split()[-1]))
       	elif	peqn in	line:
            peqn_time.append(float(line.split()[-1]))
       	elif	eeqn in	line:
            eeqn_time.append(float(line.split()[-1]))
       	elif	lagrangian in	line:
            lagrangian_time.append(float(line.split()[-1]))
    result = [execution_time,ueqn_time,peqn_time,yeqn_time,eeqn_time,lagrangian_time]
    return result

def plot_data(data):

    width = 0.5
    n = np.shape(data[0][0])[0]-1
    execution =[]
    U =[]
    pres =[]
    Y =[]
    E =[]
    lagrangian=[]
    
    for i in range(np.shape(data)[0]):
        execution.append(data[i][0][n])
        U.append(data[i][1][n])
        pres.append(data[i][2][n])
        Y.append(data[i][3][n])
        E.append(data[i][4][n])
        lagrangian.append(data[i][5][n])
    fig,ax = plt.subplots(1)
    labels = ['LoadBal + Ref','LoadBal','Standard']
    ax.bar(labels,Y,width,label='Chemistry+Transport')
    ax.bar(labels,lagrangian,width,bottom=np.array(Y),label='Lagrangian')
    ax.bar(labels,U,width,bottom=np.array(Y) + np.array(lagrangian),label = 'Momentum')
    ax.bar(labels,pres,width,bottom=np.array(Y) + np.array(lagrangian) + np.array(U),label='Pressure')
    ax.bar(labels,E,width,bottom=np.array(Y) + np.array(lagrangian) + np.array(U) + np.array(pres),label='Energy')
    rects1 = ax.bar(labels,execution,width,color='black',fill=False)
    height = rects1[0].get_height()
    ax.text(rects1[0].get_x() + rects1[0].get_width()/2., 1.1*height,r'$\chi_{su}$' + '= {:.2f}'.format(execution[2]/execution[0]),ha='center', va='bottom',fontsize=6)
    height = rects1[1].get_height()
    ax.text(rects1[1].get_x() + rects1[1].get_width()/2., 1.05*height,r'$\chi_{su}$' + '= {:.2f}'.format(execution[2]/execution[1]),ha='center', va='bottom',fontsize = 6)
    ax.legend(loc=2,frameon=False)
    ax.set_ylabel('Execution Time [s]')
    ax.ticklabel_format(axis='y',style='sci',scilimits=(2,4))
    ax.set_ylim([0,None])
    fig.tight_layout()
    fig.savefig('sprayA_bar.pdf',bbox_inches='tight',pad_inches=0)
    plt.show(block=False)
    plt.pause(2)
    plt.close()
def main():
    loadbal = get_data('../../data/3d_spray/balance_refmap/log.sprayFoam_timed')
    loadbal_noref = get_data('../../data/3d_spray/balance/log.sprayFoam_timed')
    standard = get_data('../../data/3d_spray/standard/log.sprayFoam_timed')
    plot_data([loadbal,loadbal_noref,standard])
 

main()
