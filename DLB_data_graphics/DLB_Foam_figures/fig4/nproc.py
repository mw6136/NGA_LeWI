import numpy as np
import pandas as pd
import os
import re
import matplotlib
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
    'figure.figsize': [2.64, 2.3],
    'text.usetex': True,
    'font.family': 'serif',
}

matplotlib.rcParams.update(params)

def theoretical(heavy):
    theoretical = []
    a = 1.0
    b = 0.0
    phi = 0.2
    theoretical.append(((a*heavy+b)/(((a*heavy)+(1/phi-1))/(1/phi))))

    a = 0.8
    b = 0.2
    phi = 0.25
    theoretical.append(((a*heavy+b)/(((a*heavy)+(1/phi-1))/(1/phi))))

    a = 0.4
    b = 0.6
    phi = 0.5
    theoretical.append(((a*heavy+b)/(((a*heavy)+(1/phi-1))/(1/phi))))
    
    theoretical.append(1)
    return theoretical


def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        pdb.set_trace()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),
                ha='center', va='bottom')    

def plot_bar(cases):
    fig,ax = plt.subplots()
    width = 0.1  # the width of the bars
    i = 1
    for case in cases:
        df = read_data('../../data/unit_benchmark/n_procs/' +case)
        Sr = get_Sr(df)
        h_l = get_h_l(df)
        theoretical_sr = theoretical(h_l)
        x = np.arange(len(Sr))  # the label locations
        rects1 = ax.bar(x + (width*i*1.2), Sr, width, alpha=0.7, label='$N_p$: ' + re.findall('\d+', case )[0])
        i+=1
       
    ax.axhline(y=theoretical_sr[0],xmin=0.04,xmax=0.23,color='gray',linewidth=0.8,linestyle='--')    
    ax.axhline(y=theoretical_sr[1],xmin=0.29,xmax=0.48,color='gray',linewidth=0.8,linestyle='--')    
    ax.axhline(y=theoretical_sr[2],xmin=0.53,xmax=0.72,color='gray',linewidth=0.8,linestyle='--')    
    ax.axhline(y=theoretical_sr[3],xmin=0.78,xmax=0.97,color='gray',linewidth=0.8,linestyle='--')    

    ax.set_ylabel('$\chi_{su}$')
    ax.set_xticks([0.5,1.5,2.5,3.5])
    ax.set_yticks([0,1,2,3,4])
    ax.set_xticklabels(['$C_1$','$C_2$','$C_3$','$C_4$'])
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig('speedup_ncore.pdf',bbox_inches='tight',pad_inches=0.0)
    plt.show(block=False)
    plt.pause(2)
    plt.close()

def get_Sr(df):
    grouped = df.groupby(['Function','Init. Condition'])
    Sr = []    
    x = []
    for key, item in grouped:
        func = item['Function'].unique()[0]
        cond = item['Init. Condition'].unique()[0]
        x.append(cond)
        #balancer = item['Balancer'].unique()[0]
        title = "{0}_{1}".format(func, cond)
        #print(item)
        models = item["Model"].unique()
        n_models = len(models)
        #print(models)
        for i, model in enumerate(models):
            sliced = item.loc[item["Model"] == model]
            if(i==0):
                Srtmp=np.max(sliced["Mean"])
            else:
                Srtmp = Srtmp/np.max(sliced["Mean"])
        Sr.append(Srtmp)
    return Sr
    
def get_h_l(df):

    grouped = df.groupby(['Function','Init. Condition'])
    Sr = []    
    x = []
    for key, item in grouped:
        func = item['Function'].unique()[0]
        cond = item['Init. Condition'].unique()[0]
        x.append(cond)
        title = "{0}_{1}".format(func, cond)
        models = item["Model"].unique()
        n_models = len(models)
        for i, model in enumerate(models):
            sliced = item.loc[item["Model"] == model]
            if(i==0):
                Srtmp=np.max(sliced["Mean"])/np.min(sliced["Mean"])
        Sr.append(Srtmp)

    return Sr[0]

def read_data(case):

    names = get_fnames(case)
    frames = []
    ids = []
    for name in names:      
        proc_id = int(re.search(r'\d+', name).group())
        df = pd.read_csv(case +'/' +  name)
        ddf = df.assign(Processor=proc_id)
        frames.append(ddf)
    big = pd.concat(frames, ignore_index=True)
    return big


def get_fnames(case):

    all_files = os.listdir(case+'/')
    good_files = []

    for f in all_files:
        if ("results_" in f):
            good_files.append(f)
    return good_files


def main():

    cases = ['core_40','core_80','core_160','core_320','core_640','core_1280']
    plot_bar(cases)

main()
