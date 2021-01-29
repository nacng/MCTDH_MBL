import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import csv

def scaled_wi():
    rc('text', usetex=True)
    rc('font', family='serif', serif='Computer Modern Roman')
    
    wdirlist = {'8': "./8/0.1_scaled/", 
               '12': "./12/0.1_scaled/", 
               '16': "./16/0.1_scaled/",
               '16long': "./16/0.1_scaled/longTime/",
               '24': "./24/0.1_scaled/",
               '48': "./48/0.1_scaled/",
               '96': "./96/0.1_scaled/",
               '96g0': "./96/0.1_scaled_nonint/",
               '48g0': "./48/0.1_scaled_nonint/"}
    colorList = {'8':'#4daf4a', '12':'#ff7f00', '16':'#999999', '24':'#a65628', '48':'#377eb8', '96':'#e41a1c', '96g0':'#e41a1c', '48g0':'#377eb8'}
    styleList = {'8':'solid', '12':'solid', '16':'solid', '24':'solid', '48':'solid', '96':'solid', '96g0':'dotted', '48g0':'dotted'}
    labelList = {'8':'L=8', '12':'L=12', '16':'L=16', '24':'L=24', '48':'L=48', '96':'L=96'}

    tScaleList = {'8': (0.8/np.pi)*0.01*12/8,
               '12': (0.8/np.pi)*0.01*12/12,
               '16': (0.8/np.pi)*0.01*12/16,
               '16long': (0.8/np.pi)*0.01*12/16,
               '24': (0.8/np.pi)*0.01*12/24/0.658211951,
               '48': (0.8/np.pi)*0.01*12/48/0.658211951,
               '96': (0.8/np.pi)*0.01*12/96/0.658211951}

    fig, ax = plt.subplots()
    ax.remove()
    specs = fig.add_gridspec(ncols=4, nrows=2, width_ratios=[4, 4, 4, 4], height_ratios=[3, 3])
    fig.set_figwidth(4*4)
    fig.set_figheight(6)

    yScaleList = {'8': 1/0.01, 
               '12': 1/0.01, 
               '16': 1/0.01,
               '24': 1/0.01,
               '48': 1/0.01,
               '96': 1/0.01,
               '96g0': 1/0.01,
               '48g0': 1/0.01}
    ax00 = fig.add_subplot(specs[0, 0])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"qdVar.csv", delimiter=','))
        times = c1[1:71]
        vals = c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax00.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax00.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=None)
    for n in ['16long']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"qdVar.csv", delimiter=','))
        tt = c1
        logpow = 3
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax00.errorbar(times*tScaleList['16'], vals*yScaleList['16'], stdErr*yScaleList['16'], color=colorList['16'], ecolor='k', marker='o', markeredgecolor='k', dash_capstyle='projecting', capsize=3, label=None)
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"qdVar.csv", delimiter=','))
        tt = c1
        logpow = 1.5
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax00.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax00.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=None)

    ax00.set_xscale('log')
    ax00.set_xlim(0.004, 390)
    ax00.set_ylim(0, 20.)
    ax00.set_ylabel(r'$\Delta^2_{\mathcal{Q}}/0.1^2$', fontsize=25)
    ax00.tick_params(axis='y', labelsize=17)
    ax00.tick_params(axis='x', labelsize=20)

    ax00.text(-0.15, -0.2, "\\textbf{\\framebox{(a)}}", transform=ax00.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    
######
######


    yScaleList = {'8': 8/0.01/2, 
               '12': 12/0.01/2, 
               '16': 16/0.01/2,
               '24': 24/0.01/2,
               '48': 48/0.01/2,
               '96': 96/0.01/2,
               '96g0': 96/0.01/2,
               '48g0': 48/0.01/2}
    ax01 = fig.add_subplot(specs[0, 1])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"sgQ.csv", delimiter=','))
        times = c1[1:71]
        vals = 1 - c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax01.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax01.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16long']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"sgQ.csv", delimiter=','))
        tt = c1
        logpow = 3
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax01.errorbar(times*tScaleList['16'], vals*yScaleList['16'], stdErr*yScaleList['16'], color=colorList['16'], ecolor='k', marker='o', markeredgecolor='k', dash_capstyle='projecting', capsize=3, label=None)
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"sgQ.csv", delimiter=','))
        tt = c1
        logpow = 1.5
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]
        
        ax01.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax01.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax01.set_xscale('log')
    ax01.set_xlim(0.004, 390)
    ax01.set_ylim(0, 55)
    ax01.set_ylabel(r'$\frac{L}{2}(1 - q_{EA})/0.1^2$', fontsize=25)
    ax01.tick_params(axis='y', labelsize=17)
    ax01.tick_params(axis='x', labelsize=20)

    ax01.legend(loc=9,prop={'size':15},ncol=2,frameon=True)
    ax01.text(-0.15, -0.2, "\\textbf{\\framebox{(b)}}", transform=ax01.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    
######
######


    yScaleList = {'8': 1., 
               '12': 1., 
               '16': 1.,
               '24': 1./np.log(2),
               '48': 1./np.log(2),
               '96': 1./np.log(2),
               '96g0': 1./np.log(2),
               '48g0': 1./np.log(2)}
    ax02 = fig.add_subplot(specs[0, 2])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:71]
        vals = c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax02.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax02.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16long']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"savg.csv", delimiter=','))
        tt = c1
        logpow = 3
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax02.errorbar(times*tScaleList['16'], vals*yScaleList['16'], stdErr*yScaleList['16'], color=colorList['16'], ecolor='k', marker='o', markeredgecolor='k', dash_capstyle='projecting', capsize=3, label=None)

    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"savg.csv", delimiter=','))
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

    ax02.set_xlim(0.004, 390)
    ax02.set_xscale("log")
    ax02.set_ylim(0.1,0.35)
    ax02.set_yticks([0.1, 0.2, 0.3])
    ax02.set_ylabel(r'$S_{A}$', fontsize=25)
    ax02.tick_params(axis='y', labelsize=17)
    ax02.tick_params(axis='x', labelsize=20)

    ax02.text(-0.15, -0.2, "\\textbf{\\framebox{(c)}}", transform=ax02.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')


######
######


    yScaleList = {'8': 1,
               '12': 1,
               '16': 1,
               '24': 1,
               '48': 1,
               '96': 1,
               '96g0': 1,
               '48g0': 1}
    ax03 = fig.add_subplot(specs[0, 3])
    for n in ['8', '12']:
        popData = np.genfromtxt(wdirlist[n]+"qdWF.csv", delimiter=',')
        for qdn in range(1,4):
            times = popData[1:71, 0]
            vals = 0.5*(popData[1:71, qdn] + popData[1:71, 8-qdn])
            m2L = popData[1:71, qdn+7]
            m2U = popData[1:71, 15-qdn]
            numruns = popData[0, 0]
            stdErr = 0.5*np.sqrt(m2L**2 + m2U**2)/np.sqrt(numruns)

            ax03.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
            ax03.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=None)
    for n in ['16', '24', '48', '96']:
        popData = np.genfromtxt(wdirlist[n]+"qdWF.csv", delimiter=',')
        tt = popData[0:, 0]
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        times = popData[logindx, 0]
        numruns = popData[0,0]
        for qdn in range(1,4):
            vals = 0.5*(popData[logindx, qdn] + popData[logindx, 8-qdn])
            m2L = popData[logindx, qdn+7]
            m2U = popData[logindx, 15-qdn]
            stdErr = 0.5*np.sqrt(m2L**2 + m2U**2)/np.sqrt(numruns)

            ax03.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
            ax03.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=None)

    ax03.set_xscale('log')
    ax03.set_xlim(0.004, 390)
    ax03.set_xticks([0.1,10])
    ax03.set_ylim(0, 0.04)
    ax03.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax03.yaxis.get_offset_text().set_fontsize(14)
    ax03.set_ylabel(r'$(p_n + p_{d-n})/2$', fontsize=25)
    ax03.tick_params(axis='y', labelsize=17)
    ax03.tick_params(axis='x', labelsize=20)

    ax03.text(-0.15, -0.2, "\\textbf{\\framebox{(d)}}", transform=ax03.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')

    ax03.text(110, 0.0315, "$|n-n_{\\mathrm{mid}}| = 1$", size=20, horizontalalignment='right', verticalalignment='bottom')
    ax03.text(110, 0.013, "$= 2$", size=20, horizontalalignment='right', verticalalignment='bottom')
    ax03.text(110, 0.001, "$= 3$", size=20, horizontalalignment='right', verticalalignment='bottom')
    
  
  
  
  
  
    
    idirlist = {'8': "./8/0.3_scaled/", 
               '12': "./12/0.3_scaled/", 
               '16': "./16/0.3_scaled/",
               '16long': "./16/0.3_scaled/longTime/",
               '24': "./24/0.3_scaled/",
               '48': "./48/0.3_scaled/",
               '96': "./96/0.3_scaled/",
               '96g0': "./96/0.3_scaled_nonint/",
               '48g0': "./48/0.3_scaled_nonint/"}

    tScaleList = {'8': (0.8/np.pi)*0.09*12/8,
               '12': (0.8/np.pi)*0.09*12/12,
               '16': (0.8/np.pi)*0.09*12/16,
               '16long': (0.8/np.pi)*0.09*12/16,
               '24': (0.8/np.pi)*0.09*12/24/0.658211951,
               '48': (0.8/np.pi)*0.09*12/48/0.658211951,
               '96': (0.8/np.pi)*0.09*12/96/0.658211951}
    
    
######
######


    yScaleList = {'8': 1/0.09, 
               '12': 1/0.09, 
               '16': 1/0.09,
               '16long': 1/0.09,
               '24': 1/0.09,
               '48': 1/0.09,
               '96': 1/0.09,
               '96g0': 1/0.09,
               '48g0': 1/0.09}
    ax10 = fig.add_subplot(specs[1, 0])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"qdVar.csv", delimiter=','))
        times = c1[1:71]
        vals = c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax10.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax10.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16long']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"qdVar.csv", delimiter=','))
        tt = c1
        logpow = 3
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax10.errorbar(times*tScaleList['16'], vals*yScaleList['16'], stdErr*yScaleList['16'], color=colorList['16'], ecolor='k', marker='o', markeredgecolor='k', dash_capstyle='projecting', capsize=3, label=None)
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"qdVar.csv", delimiter=','))
        tt = c1
        logpow = 1.15
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax10.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax10.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax10.set_xscale('log')
    ax10.set_xlim(0.004, 390)
    ax10.set_xlabel(r'$\gamma^2 \, t / \Omega$', fontsize=23)
    ax10.set_ylim(0, 20.0)
    ax10.set_ylabel(r'$\Delta^2_{\mathcal{Q}}/0.3^2$', fontsize=25)
    ax10.tick_params(axis='y', labelsize=17)
    ax10.tick_params(axis='x', labelsize=20)

    ax10.text(-0.15, -0.2, "\\textbf{\\framebox{(e)}}", transform=ax10.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    
######
######


    yScaleList = {'8': 8/0.09/2, 
               '12': 12/0.09/2, 
               '16': 16/0.09/2,
               '24': 24/0.09/2,
               '48': 48/0.09/2,
               '96': 96/0.09/2,
               '96g0': 96/0.09/2,
               '48g0': 48/0.09/2}
    ax11 = fig.add_subplot(specs[1, 1])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"sgQ.csv", delimiter=','))
        times = c1[1:71]
        vals = 1 - c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax11.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax11.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16long']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"sgQ.csv", delimiter=','))
        tt = c1
        logpow = 3
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1- c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax11.errorbar(times*tScaleList['16'], vals*yScaleList['16'], stdErr*yScaleList['16'], color=colorList['16'], ecolor='k', marker='o', markeredgecolor='k', dash_capstyle='projecting', capsize=3, label=None)
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"sgQ.csv", delimiter=','))
        tt = c1
        logpow = 1.15
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax11.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax11.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax11.set_xscale('log')
    ax11.set_xlim(0.004, 390)
    ax11.set_xlabel(r'$\gamma^2 \, t / \Omega$', fontsize=23)
    ax11.set_ylim(0, 55)
    ax11.set_ylabel(r'$\frac{L}{2}(1 - q_{EA})/0.3^2$', fontsize=25)
    ax11.tick_params(axis='y', labelsize=17)
    ax11.tick_params(axis='x', labelsize=20)

    ax11.text(-0.15, -0.2, "\\textbf{\\framebox{(f)}}", transform=ax11.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    
######
######


    yScaleList = {'8': 1., 
               '12': 1., 
               '16': 1.,
               '24': 1./np.log(2),
               '48': 1./np.log(2),
               '96': 1./np.log(2),
               '96g0': 1./np.log(2),
               '48g0': 1./np.log(2)}
    ax12 = fig.add_subplot(specs[1, 2])
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:71]
        vals = c2[1:71]
        m2 = c3[1:71]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax12.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax12.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['16long']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"savg.csv", delimiter=','))
        tt = c1
        logpow = 3
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax12.errorbar(times*tScaleList['16'], vals*yScaleList['16'], stdErr*yScaleList['16'], color=colorList['16'], ecolor='k', marker='o', markeredgecolor='k', dash_capstyle='projecting', capsize=3, label=None)
    for n in ['16', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"savg.csv", delimiter=','))
        tt = c1
        logpow = 1.15
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax12.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax12.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    tvar = np.logspace(-1, 3, base=10, num=20)
    print(tvar)
    yvar = [np.log(t + 0.6) + 1.7 for t in tvar]
    ax12.plot(tvar, yvar, linewidth=1.5, marker='x', color = 'k', linestyle="dotted", label=r'$\log(0.6 + \frac{\gamma^2}{\Omega} t)+1.7$')
    ax12.legend(loc=4,prop={'size':15},frameon=False)

    ax12.set_xlim(0.004, 390)
    ax12.set_xlabel(r'$\gamma^2 \, t / \Omega$', fontsize=23)
    ax12.set_ylim(0,5.65)
    ax12.set_xscale('log')
    ax12.set_yticks([0, 1, 2, 3, 4, 5])
    ax12.set_ylabel(r'$S_{A}$', fontsize=25)
    ax12.tick_params(axis='y', labelsize=17)
    ax12.tick_params(axis='x', labelsize=20)

    ax12.text(-0.15, -0.2, "\\textbf{\\framebox{(g)}}", transform=ax12.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    
######
######


    yScaleList = {'8': 1,
               '12': 1,
               '16': 1,
               '24': 1,
               '48': 1,
               '96': 1,
               '96g0': 1,
               '48g0': 1}
    ax13 = fig.add_subplot(specs[1, 3])
    for n in ['8', '12']:
        popData = np.genfromtxt(idirlist[n]+"qdWF.csv", delimiter=',')
        for qdn in range(1,4):
            times = popData[1:71, 0]
            vals = 0.5*(popData[1:71, qdn] + popData[1:71, 8-qdn])
            m2L = popData[1:71, qdn+7]
            m2U = popData[1:71, 15-qdn]
            numruns = popData[0,0]
            stdErr = 0.5*np.sqrt(m2L**2 + m2U**2)/np.sqrt(numruns)

            ax13.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
            ax13.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])
    for n in ['16', '24', '48', '96']:
        popData = np.genfromtxt(idirlist[n]+"qdWF.csv", delimiter=',')
        tt = popData[0:, 0]
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(tt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        times = popData[logindx, 0]
        numruns = popData[0,0]
        for qdn in range(1,4):
            vals = 0.5*(popData[logindx, qdn] + popData[logindx, 8-qdn])
            m2L = popData[logindx, qdn+7]
            m2U = popData[logindx, 15-qdn]
            stdErr = 0.5*np.sqrt(m2L**2 + m2U**2)/np.sqrt(numruns)

            ax13.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
            ax13.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=labelList[n])

    ax13.set_xscale('log')
    ax13.set_xlim(0.004, 390)
    ax13.set_xticks([0.1,10])
    ax13.set_ylim(0, 0.2)
    ax13.set_xlabel(r'$\gamma^2 \, t / \Omega$', fontsize=23)
    ax13.set_ylim(0, 0.23)
    ax13.set_yticks([0,0.1,0.2])
    ax13.set_ylabel(r'$(p_n + p_{d-n})/2$', fontsize=25)
    ax13.tick_params(axis='y', labelsize=17)
    ax13.tick_params(axis='x', labelsize=20)

    ax13.text(-0.15, -0.2, "\\textbf{\\framebox{(h)}}", transform=ax13.transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    ax13.text(110, 0.2, "$|n-n_{\\mathrm{mid}}| = 1$", size=20, horizontalalignment='right', verticalalignment='center')
    ax13.text(110, 0.1, "$= 2$", size=20, horizontalalignment='right', verticalalignment='center')
    ax13.text(110, 0.03, "$= 3$", size=20, horizontalalignment='right', verticalalignment='center')

    
    fig.tight_layout()
    fig.savefig('dynamics_scaled.pdf')
