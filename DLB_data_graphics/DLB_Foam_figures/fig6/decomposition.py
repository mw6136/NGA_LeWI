import matplotlib.pyplot as plt
import numpy as np
import helper as vt
import matplotlib
import rand_cmap as rnd
params = {
    'image.origin': 'lower',
    'image.interpolation': 'nearest',
    'image.cmap': 'gray',
    'axes.grid': False,
    'savefig.dpi': 1000,  # to adjust notebook inline plot size
    'axes.labelsize': 9, # fontsize for x and y labels (was 10)
    'axes.titlesize': 9,
    'font.size': 9, # was 10
    'legend.fontsize': 9, # was 10
    'xtick.labelsize': 9,
    'text.usetex': True,
    'ytick.labelsize': 9,
    'figure.figsize': [5.64, 5],
    'font.family': 'serif',
}
matplotlib.rcParams.update(params)

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return matplotlib.colors.ListedColormap(color_list, name=cmap_name)


cases = ['simple_x/loadbal_32/','simple_y/loadbal_32/','scotch/loadbal_32/']
labels = ['Simple-X','Simple-Y','Scotch']

fig,ax = plt.subplots(nrows=1,ncols=np.size(cases))

for idx,case in enumerate(cases):

    path ='../../data/shearlayer/'+case+'postProcessing/decomposition/0/cellDist_zNormal.raw'

    x_coords = []
    y_coords = []
    z_coords = []
    cellDist = []
    with open(path) as f:
        next(f)
        next(f)
        for line in f:
            data = line.split()
            x_coords.append(float(data[0]))
       	    y_coords.append(float(data[1]))
       	    z_coords.append(float(data[2]))
       	    cellDist.append(int(data[3]))

    xi,yi = vt.generate_target_grid(x_coords,y_coords,20e-6,20e-6)
    cellDist_interp = vt.interpolate_grid(x_coords,y_coords,cellDist,xi,yi)
    N = np.max(cellDist) - np.min(cellDist)
    im = ax[idx].contourf(xi,yi,cellDist_interp,np.linspace(np.min(cellDist),np.max(cellDist),np.max(cellDist) - np.min(cellDist)),cmap=rnd.rand_cmap(N, type='soft'),extend='both')
    ax[idx].set_xticks([])
    ax[idx].set_yticks([])
    ax[idx].set_aspect(1)
    ax[idx].set_title(labels[idx])

cbar = fig.colorbar(im, ax=ax[:],orientation='horizontal',fraction=0.02, pad=0.015,ticks=[np.min(cellDist),np.max(cellDist)])
cbar.ax.set_xlabel('Processor ID')
fig.savefig('decomposition.png',bbox_inches='tight',pad_inches=0.0)   
plt.show(block=False)
plt.pause(2)
plt.close()