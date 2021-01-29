import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import csv

def short_time_wi():
    rc('text', usetex=True)
    rc('font', family='serif', serif='Computer Modern Roman')
    
    wdirlist = {'8': "./8/0.1_scaled/", 
               '12': "./12/0.1_scaled/", 
               '16': "./16/0.1_scaled/",
               '24': "./24/0.1_scaled/",
               '48': "./48/0.1_scaled/",
               '96': "./96/0.1_scaled/",
               '96g0': "./96/0.1_scaled_nonint/",
               '48g0': "./48/0.1_scaled_nonint/"}
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

    tt, qd_var_avg = np.transpose(np.genfromtxt("./pertThQdVar.csv", delimiter=","))

    #style/label parameters
    fig, ax = plt.subplots()
    ax.remove()
    specs = fig.add_gridspec(ncols=3, nrows=2, width_ratios=[5, 5, 5], height_ratios=[3, 3])
    fig.set_figwidth(15)
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
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"qdVar.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax00.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax00.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=None)
    for n in ['48g0', '96g0', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"qdVar.csv", delimiter=','))
        ttt = c1
        logpow = 1.15
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax00.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax00.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n], label=None)

    ax00.plot(tt*tScaleList['12'], qd_var_avg / 0.01, label="Perturbation theory", color = 'k', linestyle="dotted")
    legend_elements = [
                       Line2D([0],[0], color=colorList['12'], ls='-', label=labelList['12']),
                       Line2D([0],[0], color=colorList['24'], ls='-', label=labelList['24']),
                       Line2D([0],[0], color=colorList['48'], ls='-', label=labelList['48']),
                       Line2D([0],[0], color=colorList['96'], ls='-', label=labelList['96']),
                      ]
    ax00.legend(handles=legend_elements,loc=2,prop={'size':15},ncol=2,frameon=False)

    ax00.set_xscale('log')
    ax00.set_xlim(0.1, 160)
    ax00.set_ylabel(r'$\Delta^2_{\mathcal{Q}}/0.1^2$', fontsize=25)
    ax00.tick_params(axis='y', labelsize=17)
    ax00.tick_params(axis='x', labelsize=20)

    ax00.text(-0.18, -0.2, "\\textbf{\\framebox{(a)}}", transform=ax00.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    

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
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"sgQ.csv", delimiter=','))
        times = c1[1:41]
        vals = 1 - c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax01.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax01.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['48g0', '96g0', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"sgQ.csv", delimiter=','))
        ttt = c1
        logpow = 1.05
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax01.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax01.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    ax01.plot(tt*tScaleList['12'], qd_var_avg / 0.01, label="Perturbation theory", color = 'k', linestyle="dotted")

    ax01.set_xlim(0.1, 160)
    ax01.set_xscale("log")
    ax01.set_ylim(0, 40)
    ax01.set_ylabel(r'$\frac{L}{2}(1 - q_{EA})/0.1^2$', fontsize=25)
    ax01.tick_params(axis='y', labelsize=17)
    ax01.tick_params(axis='x', labelsize=20)

    ax01.text(-0.18, -0.2, "\\textbf{\\framebox{(b)}}", transform=ax01.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    

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
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax02.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax02.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['48g0', '96g0', '24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(wdirlist[n]+"savg.csv", delimiter=','))
        ttt = c1
        logpow = 1.05
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax02.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax02.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    ax02.set_xlim(0.1, 160)
    ax02.set_xscale("log")
    ax02.set_ylabel(r'$S_{A}$', fontsize=25)
    ax02.tick_params(axis='y', labelsize=17)
    ax02.tick_params(axis='x', labelsize=20)

    ax02.text(-0.18, -0.2, "\\textbf{\\framebox{(c)}}", transform=ax02.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')










    idirlist = {'8': "./8/0.3_scaled/", 
               '12': "./12/0.3_scaled/", 
               '16': "./16/0.3_scaled/",
               '24': "./24/0.3_scaled/",
               '48': "./48/0.3_scaled/",
               '96': "./96/0.3_scaled/",
               '96g0': "./96/0.3_scaled_nonint/",
               '48g0': "./48/0.3_scaled_nonint/"}

    yScaleList = {'8': 1/0.09,
               '12': 1/0.09,
               '16': 1/0.09,
               '24': 1/0.09,
               '48': 1/0.09,
               '96': 1/0.09,
               '96g0': 1/0.09,
               '48g0': 1/0.09}
    ax10 = fig.add_subplot(specs[1, 0])
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"qdVar.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax10.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax10.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"qdVar.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax10.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax10.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['96g0', '48g0']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"qdVar.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax10.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax10.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    ax10.plot(tt*tScaleList['12'], qd_var_avg / 0.01, label="Perturbation theory", color = 'k', linestyle="dotted")

    ax10.set_xscale('log')
    ax10.set_xlim(0.1, 160)
    ax10.set_xlabel('$\Omega \, t/(2 \pi)$', fontsize=23)
    ax10.set_ylabel(r'$\Delta^2_{\mathcal{Q}}/0.3^2$', fontsize=25)
    ax10.tick_params(axis='y', labelsize=17)
    ax10.tick_params(axis='x', labelsize=20)

    ax10.legend(handles=[Line2D([0],[0],color='k',ls="--",label="Noninteracting spins ($g=0$)"),Line2D([0],[0],color='k',ls='dotted',label="Perturbation theory")],loc=2,prop={'size':15},ncol=1,frameon=False)

    ax10.text(-0.18, -0.2, "\\textbf{\\framebox{(d)}}", transform=ax10.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    

######
######


    yScaleList = {'8': 0.5*8/0.09,
               '12': 0.5*12/0.09,
               '16': 0.5*16/0.09,
               '24': 0.5*24/0.09,
               '48': 0.5*48/0.09,
               '96': 0.5*96/0.09,
               '96g0': 0.5*96/0.09,
               '48g0': 0.5*48/0.09}
    ax11 = fig.add_subplot(specs[1, 1])
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"sgQ.csv", delimiter=','))
        times = c1[1:41]
        vals = 1 - c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = np.sqrt(m2)/np.sqrt(numruns*(numruns-1))

        ax11.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax11.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"sgQ.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax11.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax11.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['96g0', '48g0']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"sgQ.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = 1 - c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]
        
        ax11.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax11.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    ax11.plot(tt*tScaleList['12'], qd_var_avg / 0.01, label="Perturbation theory", color = 'k', linestyle="dotted")

    ax11.set_xlim(0.1, 160)
    ax11.set_xscale("log")
    ax11.set_xlabel('$\Omega \, t/(2 \pi)$', fontsize=23)
    ax11.set_ylim(0, 40)
    ax11.set_ylabel(r'$\frac{L}{2} (1 - q_{EA}) /0.3^2$', fontsize=25)
    ax11.tick_params(axis='y', labelsize=17)
    ax11.tick_params(axis='x', labelsize=20)

    ax11.text(-0.18, -0.2, "\\textbf{\\framebox{(e)}}", transform=ax11.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    

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
    for n in ['12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax12.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.15)
        ax12.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['24', '48', '96']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"savg.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax12.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax12.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['96g0', '48g0']:
        c1, c2, c3 = np.transpose(np.genfromtxt(idirlist[n]+"savg.csv", delimiter=','))
        ttt = c1
        logpow = 1.1
        logindx = np.unique(np.floor(np.logspace(1, np.log(len(ttt) - 1)/np.log(logpow), base=logpow, endpoint=True)).astype('int'))
        vals = c2[logindx]
        m2 = c3[logindx]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)
        times = c1[logindx]

        ax12.fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax12.plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])

    ax12.set_xlim(0.1, 160)
    ax12.set_xlabel('$\Omega \, t/(2 \pi)$', fontsize=23)
    ax12.set_ylim(0,4.0)
    ax12.set_xscale('log')
    ax12.set_yticks([0, 1, 2, 3, 4])
    ax12.set_ylabel(r'$S_{A}$', fontsize=25)
    ax12.tick_params(axis='y', labelsize=17)
    ax12.tick_params(axis='x', labelsize=20)

    ax12.text(-0.18, -0.2, "\\textbf{\\framebox{(f)}}", transform=ax12.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='left')
    
    fig.tight_layout()
    fig.savefig('dynamics_wi.pdf')
