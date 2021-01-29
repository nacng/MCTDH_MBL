import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib.lines import Line2D
import matplotlib.patches as patches
import csv

def dynamics_strong():
    rc('text', usetex=True)
    rc('font', family='serif', serif='Computer Modern Roman')
    
    dirlist = {'8': "./8/0.7_scaled/", 
               '12': "./12/0.7_scaled/", 
               '16': "./16/0.7_scaled/",
               '20': "./20/0.7_scaled/", 
               '24cheb': "./24_Chebyshev/0.7_scaled/", 
               '24': "./24/0.7_scaled/",
               '48': "./48/0.7_scaled/",
               '96': "./96/0.7_scaled/"}
    colorList = {'8':'#4daf4a', '12':'#ff7f00', '16':'#999999', '20':"purple", '24cheb': '#a65628', '24':'#a65628', '48':'#377eb8', '96':'#e41a1c', '96g0':'#e41a1c', '48g0':'#377eb8', '80': '#e41a1c'}
    styleList = {'8':'dashdot', '12':'dashdot', '16':'dashdot', '20':"dashdot", '24cheb': 'dashdot', '24':'solid', '48':'solid', '96':'solid', '80': 'solid', '96g0':'dotted', '48g0':'dotted'}
    labelList = {'8':'L=8', '12':'L=12', '16':'L=16', '20':'L=20', '24cheb':'L=24', '24':'L=24', '48':'L=48', '96':'L=96', '80':'L=80'}

    tScaleList = {'8': (1/1.6), 
               '12': (1/1.6),
               '16': (1/1.6),
               '20': (1/1.6),
               '24cheb': (1/1.6),
               '24': (1/1.6)/0.658211951,
               '48': (1/1.6)/0.658211951,
               '80': (1/1.6)/0.658211951}

    ent_lim_24 = 112.0
    ent_lim_48 = 192.0
    ent_lim_80 = 160.0
    entanglement_limit = {'24':np.log(ent_lim_24), '48':np.log(ent_lim_48), '80':np.log(ent_lim_80)}


    t_24 = np.loadtxt('strong_not_converged/24/time.dat')
    t_48 = np.loadtxt('strong_not_converged/48/time.dat')
    t_80 = np.loadtxt('strong_not_converged/80/time.dat')

    q_24 = np.loadtxt('strong_not_converged/24/spin_glass_exp.dat')
    q_48 = np.loadtxt('strong_not_converged/48/spin_glass_exp.dat')
    q_80 = np.loadtxt('strong_not_converged/80/spin_glass_exp.dat')

    qd_var_24 = np.loadtxt('strong_not_converged/24/qudit_var.dat')
    qd_var_48 = np.loadtxt('strong_not_converged/48/qudit_var.dat')
    qd_var_80 = np.loadtxt('strong_not_converged/80/qudit_var.dat')

    ent_24 = np.loadtxt('strong_not_converged/24/entanglement_A.dat')
    ent_48 = np.loadtxt('strong_not_converged/48/entanglement_A.dat')
    ent_80 = np.loadtxt('strong_not_converged/80/entanglement_A.dat')

    data = {'24': (t_24, qd_var_24, q_24, ent_24),
            '48': (t_48, qd_var_48, q_48, ent_48),
            '80': (t_80, qd_var_80, q_80, ent_80)}

    fig, ax = plt.subplots(ncols=2)
    fig.set_figheight(4.5)
    fig.set_figwidth(8)

    yScaleList = {'8': 1., 
               '12': 1., 
               '16': 1.,
               '20': 1.,
               '24cheb': 1.,
               '24': 1.,
               '48': 1.,
               '80': 1.,
               '96g0': 1.,
               '48g0': 1.}
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"qdVar.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax[0].fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax[0].plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['16', '20', '24cheb']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"qdVar.csv", delimiter=','))
        times = c1[1:]
        vals = c2[1:]
        m2 = c3[1:]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax[0].fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax[0].plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])


    ax[0].set_xlim(0, 30)
    ax[0].set_xlabel('$\Omega t/(2\pi)$', fontsize=25)
    ax[0].set_ylim(1.8, 4.2)
    ax[0].set_yticks([0, 1.0, 2.0, 3.0, 4.0])
    ax[0].set_ylabel(r'$\Delta^2_{\mathcal{Q}}$', fontsize=25)
    ax[0].tick_params(axis='y', labelsize=20)
    ax[0].tick_params(axis='x', labelsize=23)
    legend_elements = [
                       Line2D([0],[0], color=colorList['8'], ls='-', label=labelList['8']),
                       Line2D([0],[0], color=colorList['12'], ls='-', label=labelList['12']),
                       Line2D([0],[0], color=colorList['16'], ls='-', label=labelList['16']),
                       Line2D([0],[0], color=colorList['20'], ls='-', label=labelList['20']),
                       Line2D([0],[0], color=colorList['24'], ls='-', label=labelList['24']),
                       Line2D([0],[0], color=colorList['48'], ls='-', label=labelList['48']),
                       Line2D([0],[0], color=colorList['80'], ls='-', label=labelList['80'])
                      ]
    ax[0].legend(handles=legend_elements,loc=9,prop={'size':15},ncol=2,frameon=False)
    ax[0].text(-0.15, -0.2, "\\textbf{\\framebox{(a)}}", transform=ax[0].transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    

######
######


    yScaleList = {'8': 1., 
               '12': 1., 
               '16': 1.,
               '20': 1.,
               '24cheb': 1.,
               '24': 1./np.log(2),
               '48': 1./np.log(2),
               '80': 1./np.log(2)}
    for n in ['8', '12']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:41]
        vals = c2[1:41]
        m2 = c3[1:41]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax[1].fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax[1].plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['16']:
        c1, c2, c3 = np.transpose(np.genfromtxt(dirlist[n]+"savg.csv", delimiter=','))
        times = c1[1:]
        vals = c2[1:]
        m2 = c3[1:]
        numruns = c1[0]
        stdErr = m2/np.sqrt(numruns)

        ax[1].fill_between(times*tScaleList[n], (vals+stdErr)*yScaleList[n], (vals-stdErr)*yScaleList[n], facecolor=colorList[n], edgecolor=colorList[n], alpha=0.2)
        ax[1].plot(times*tScaleList[n], vals*yScaleList[n], linewidth=1.5, color=colorList[n], linestyle=styleList[n])
    for n in ['24', '48', '80']:
        ax[1].plot(data[n][0], data[n][3]*yScaleList[n], label=labelList[n], color = colorList[n])
        ax[1].hlines(entanglement_limit[n]*yScaleList[n], 0, 50, color = colorList[n], linestyle = 'dotted')
#    ax[1].legend(loc=4,prop={'size':20})
    ax[1].set_xlim(0, 40)
    ax[1].set_xlabel(r'$\Omega t/(2\pi)$', fontsize=25)
    ax[1].set_ylim(0, 8.5)
    ax[1].set_yticks([0, 2, 4, 6, 8])
    ax[1].set_ylabel(r'$S_A$', fontsize=25)
    ax[1].tick_params(axis='y', labelsize=20)
    ax[1].tick_params(axis='x', labelsize=23)
    ax[1].text(-0.15, -0.2, "\\textbf{\\framebox{(b)}}", transform=ax[1].transAxes, font='Arial Black', fontsize=16, fontweight='bold', va='bottom', ha='left')
    

    fig.tight_layout()
    fig.savefig('dynamics_strong.pdf')
