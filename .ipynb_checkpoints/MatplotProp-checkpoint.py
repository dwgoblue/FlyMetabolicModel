##########################Packages############################

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import seaborn as sns

##############################################################
def CanvasStyle(ax, square=False):
    ax.patch.set_facecolor('white')
    ax.grid('off')
    plt.tight_layout()
    ax.grid(False)
    ax.yaxis.set_tick_params(width=2)
    ax.xaxis.set_tick_params(width=2)
    
    if square==False:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.axhline(linewidth=3, y=ax.get_ylim()[0], color='k')
        ax.axvline(linewidth=3, x=ax.get_xlim()[0], color='k')

    else:
        ax.axhline(linewidth=3, y=ax.get_ylim()[0], color='k')
        ax.axvline(linewidth=3, x=ax.get_xlim()[0], color='k')
        ax.axhline(linewidth=3, y=ax.get_ylim()[1], color='k')
        ax.axvline(linewidth=3, x=ax.get_xlim()[1], color='k')

    return ax

def PltProps():
    plt.rcParams['font.weight'] = 'normal'#'bold'
    plt.rcParams['font.size'] = 20
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['axes.labelsize'] = 20
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams["lines.linewidth"] = 3
    plt.rcParams["axes.titlesize"] = 42
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    errbarattr = {'lw':5,
                       'capthick':3,
                       'capsize':10,
                       'ecolor':'black',
                       'markersize':15,
                       'marker':'o'
                      }

def SeabornProp():
    style = {'axes.facecolor': 'white',
 'axes.edgecolor': '.8',
 'axes.grid': False,
 'axes.axisbelow': True,
 'axes.labelcolor': '.15',
 'figure.facecolor': 'white',
 'grid.color': '.8',
 'grid.linestyle': '-',
 'text.color': '.15',
 'xtick.color': '.15',
 'ytick.color': '.15',
 'xtick.direction': 'out',
 'ytick.direction': 'out',
 'lines.solid_capstyle': 'round',
 'patch.edgecolor': 'w',
 'image.cmap': 'rocket',
 'font.family': ['sans-serif'],
 'font.sans-serif': ['Arial',
  'DejaVu Sans',
  'Liberation Sans',
  'Bitstream Vera Sans',
  'sans-serif'],
 'patch.force_edgecolor': True,
 'xtick.bottom': False,
 'xtick.top': False,
 'ytick.left': False,
 'ytick.right': False,
 'axes.spines.left': False,
 'axes.spines.bottom': True,
 'axes.spines.right': False,
 'axes.spines.top': True,}
    return style
##############################################################