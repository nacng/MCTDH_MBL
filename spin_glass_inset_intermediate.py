import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import csv

def q_inset_int():
    rc('text', usetex=True)
    rc('font', family='serif', serif='Computer Modern Roman')
    
    dirlist = {'8': "./8/0.3_scaled/", 
               '12': "./12/0.3_scaled/", 
               '16': "./16/0.3_scaled/",
               '24': "./24/0.3_scaled/",
               '48': "./48/0.3_scaled/",
               '96': "./96/0.3_scaled/",
               '96g0': "./96/0.3_scaled_nonint/",
               '48g0': "./48/0.3_scaled_nonint/"}
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
    fig.set_figwidth(5)
    fig.set_figheight(3)
               
    yScaleList = {'8': 1, 
               '12': 1, 
               '16': 1,
               '24': 1,
               '48': 1,
               '96': 1,
               '96g0': 1,
               '48g0': 1}
    ax = fig.add_subplot(specs[0, 0])
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"sgQ.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"sgQ.csv", delimiter=','))
        ttt = c1
        logpow = 1.5
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow) , base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

        ax.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    ax.set_xlim(0.1, 160)
    ax.set_xscale('log')
    ax.set_xlabel("$\\Omega t/(2\\pi)$", fontsize=30)
    ax.xaxis.set_label_coords(0.35, -0.3)
    ax.set_yticks([0.6, 0.8, 1])
    ax.set_ylabel(r'$q_{EA}$', fontsize=35)
    ax.tick_params(axis='y', labelsize=30)
    ax.tick_params(axis='x', labelsize=30)
   
    fig.tight_layout()
    fig.savefig('qInset_int.pdf', transparent=True)
