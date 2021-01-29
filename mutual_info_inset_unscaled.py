import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import csv

def mi_inset_unscaled():
    rc('text', usetex=True)
    rc('font', family='serif', serif='Computer Modern Roman')
    
    dirlist = {'8': "./8/0.10607_unscaled/", 
               '12': "./12/0.10607_unscaled/", 
               '16': "./16/0.10607_unscaled/",
               '24': "./24/0.10607_unscaled/",
               '48': "./48/0.10607_unscaled/",
               '96': "./96/0.3_scaled/",
               '96g0': "./96/0.3_scaled_nonint/"}
    colorList = {'8':'#4daf4a', '12':'#ff7f00', '16':'#999999', '24':'#a65628', '48':'#377eb8', '96':'#e41a1c', '96g0':'#e41a1c', '48g0':'#377eb8'}
    styleList = {'8':'solid', '12':'solid', '16':'solid', '24':'solid', '48':'solid', '96':'solid', '96g0':'dotted', '48g0':'dotted'}
    labelList = {'8':'L=8', '12':'L=12', '16':'L=16', '24':'L=24', '48':'L=48', '96':'L=96'}

    tScaleList = {'8': (1/1.6),
               '12': (1/1.6),
               '16': (1/1.6),
               '24': (1/1.6)/0.658211951,
               '48': (1/1.6)/0.658211951,
               '96': (1/1.6)/0.658211951,
               '96g0': (1/1.6)/0.658211951,
               '48g0': (1/1.6)/0.658211951}

    fig, ax = plt.subplots()
    ax.remove()
    specs = fig.add_gridspec(ncols=1, nrows=1, width_ratios=[1], height_ratios=[1])
    fig.set_figwidth(3)
    fig.set_figheight(3)

    yScaleList = {'8': 1./8., 
               '12': 1./12., 
               '16': 1./16.,
               '24': 1./np.log(2)/24.,
               '48': 1./np.log(2)/48.,
               '96': 1./np.log(2)/96.,
               '96g0': 1./np.log(2)/96.}
    ax = fig.add_subplot(specs[0, 0])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"MIavg.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"MIavg.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax.set_xlim(0.1, 160)
    ax.set_xscale("log")
#    ax.set_xlabel("$\\Omega t/(2\\pi)$", fontsize=30)
    ax.yaxis.set_label_coords(-0.1, 0.5)
    ax.set_ylim(0, 0.01)
    ax.set_yticks([0, 0.01])
    ax.set_ylabel(r'$\mathrm{MI}/L$', fontsize=35)
    ax.tick_params(axis='y', labelsize=30)
    ax.tick_params(axis='x', labelsize=30)

    fig.tight_layout()
    fig.savefig('MIinset_unscaled.pdf', transparent=True)
