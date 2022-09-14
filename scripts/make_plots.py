"""make various plots for FRBs
Jessie Thwaites 
September 2022
"""

import numpy as np
import matplotlib as mpl
mpl.use('agg')

import csky as cy
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import healpy as hp
import histlite as hl
import meander
import matplotlib.pyplot as plt

import sys
frb_scripts_path='/home/jthwaites/FRB/scripts'
if frb_scripts_path not in sys.path:
    sys.path.append(frb_scripts_path)
import plot_functions as plotf

plotf.sens_tw_plot('FRB110523',[1000.,10000.])
