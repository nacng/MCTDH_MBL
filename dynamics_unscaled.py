import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import csv

def fixed_coupling():
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
    specs = fig.add_gridspec(ncols=2, nrows=1, width_ratios=[5, 5])
    fig.set_figwidth(7)
    fig.set_figheight(5)

    yScaleList = {'8': 1,
               '12': 1,
               '16': 1,
               '24': 1,
               '48': 1,
               '96': 1,
               '96g0': 1,
               '48g0': 1}
    ax01 = fig.add_subplot(specs[0, 0])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"sgQ.csv", delimiter=','))
        times = c1[1:71]
        vals = 1 - c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax01.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax01.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"sgQ.csv", delimiter=','))
        tt = c1
        logpow = 1.15
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax01.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax01.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax01.set_xlim(0.1, 1000)
    ax01.set_xscale('log')
    ax01.set_xlabel(r'$\Omega t/(2\pi)$', fontsize=25)
    ax01.set_ylim(0, 0.05)
    ax01.set_yticks([0, 0.01, 0.02, 0.03, 0.04, 0.05])
    ax01.set_ylabel(r'$1 - q_{EA}$', fontsize=25)
    ax01.tick_params(axis='y', labelsize=17)
    ax01.tick_params(axis='x', labelsize=20)
    ax01.legend(loc=2,prop={'size':14},ncol=1,frameon=False)

    ax01.text(-0.3, -0.2, "\\textbf{\\framebox{(a)}}", transform=ax01.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    
######
######


    yScaleList = {'8': 1./8, 
               '12': 1./12, 
               '16': 1./16,
               '24': 1./np.log(2)/24,
               '48': 1./np.log(2)/48,
               '96': 1./np.log(2)/96,
               '96g0': 1./np.log(2)/96,
               '48g0': 1./np.log(2)/48}
    ax02 = fig.add_subplot(specs[0, 1])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:71]
        vals = c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax02.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax02.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"savg.csv", delimiter=','))
        tt = c1
        logpow = 1.15
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax02.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax02.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax02.set_xlim(0.1, 1000)
    ax02.set_xscale('log')
    ax02.set_xlabel(r'$\Omega t/(2\pi)$', fontsize=25)
    ax02.set_ylim(-0.01,0.035)
    ax02.set_yticks([0, 0.01, 0.02, 0.03])
    ax02.set_ylabel(r'$S_{A}/L$', fontsize=25)
    ax02.tick_params(axis='y', labelsize=17)
    ax02.tick_params(axis='x', labelsize=20)
    ax02.text(-0.3, -0.2, "\\textbf{\\framebox{(b)}}", transform=ax02.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    fig.tight_layout()
    fig.savefig('fixed_coupling.pdf')
