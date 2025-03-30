#-----------------------------------------------------------------------#
# castl.mcphot v0.6.2
# By Hunter Brooks, at NAU, Flagstaff: Mar. 26, 2025
#
# Purpose: Perform MCMC calculation on model spectral photometry
#
# There were 2 developers, 
# myself and the God, 
# and I have forgotten how it works.
#-----------------------------------------------------------------------#

# ------ Import Markov-Chain Monte-Carlo Packages ------ # 
import emcee
# ------------------------------------------------------- # 

# ------ Import File Processing Packages ------ #
from IPython.display import clear_output
from tqdm import tqdm
import h5py
# --------------------------------------------- # 

# ------ Import File Loading Packages ------ #
from astropy.table import Table
from astropy.io import fits
import pandas as pd
# ------------------------------------------ #

# ------ Import Plotting Packages ------ #
from IPython.display import display, Math
import matplotlib.pyplot as plt
import corner
# -------------------------------------- #

# ------ Import Math Packages ------ #
from scipy.interpolate import LinearNDInterpolator
from sklearn.preprocessing import MinMaxScaler
from scipy.interpolate import RBFInterpolator
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import astropy.units as u
import numpy as np
# ---------------------------------- #

# ------ Ignore all warnings ------ #
import warnings
warnings.filterwarnings('ignore')
# --------------------------------- #

# --------------------------------- #
def mcphot(input_file, output_file, model_directory, model_parm, 
          grid_scale = 10, unit_wave=[u.um, u.um], unit_flux=[(u.erg / (u.cm**2 * u.s * u.um)), (u.erg / (u.cm**2 * u.s * u.um))], 
          walkers=15, steps=1000, 
          rv_fit=False, monitor=False, save_output=True): 
    
    # JUST MAKE LIKE MCSPEC PRETTY MUCH
    
    return 0
    # ---------------------------------------------------- #
# --------------------------------- #

# --------------------------------- #
def obcolor(TEST): 

    # READ IN PHOTOMETRIC COLORS

    return np.nan

# --------------------------------- #
def intercolor(TEST):
    
    # READ IN MODEL AND CREATE PHOTOMETRY
    
    # COMPARE OBSERVED PHOTOMETRY AND MODEL PHOTOMETRY FOR BEST POINT
    
    # CREATE GRID AROUND BEST POINT LIKE DONE IN MCSPEC
    
    # CREATE INTERPOLATOR WITH REDUCED GRID
    
    return 0
# --------------------------------- #

# --------------------------------- #
def statmc(TEST):
    
    # CALCULATE NEW INTERPOLATED COLOR
    
    # CALCULATE CHI SQUARE FOR THE COLOR AND TAKE A NORMALIZED AVERAGE

    return 0
# --------------------------------- #

# --------------------------------- #
def prior(parm, parm_bound):
    if np.any(parm < np.array([low for low, high in parm_bound])) or np.any(parm > np.array([high for low, high in parm_bound])):
        return -np.inf
    return 0
# --------------------------------- #

# --------------------------------- #
def log_posterior(parm, observed_wave, observed_flux, unc, interpolator, scaler, parm_bound, rv_fit):
    lp = prior(parm, parm_bound)
    if not np.isfinite(lp):
        return -np.inf
    
    stat = statmc(observed_wave, observed_flux, unc, interpolator, scaler, parm, rv_fit)  
    return lp + stat if stat != 0 else -np.inf
# --------------------------------- #

# --------------------------------- #
def photmc(TEST):
    
    # SET UP INITIAL POSITIONS
    
    # SET UP SAMPLER
    
    # START RUNNING AND A MONITOR WITH GAMMA STEPA
    return 0
# --------------------------------- #

# --------------------------------- #
def mcbest(TEST): 
    
    
    # PRINT OUT BEST GRID POINT
    
    # CALCULATE CHI SQUARE AND CALCULATE
    
    return 0
# --------------------------------- #

# --------------------------------- #
def mcplot(TEST): 
    
    # CORNER PLOT
    
    # STEPPA PLOT
    
    return 0
# --------------------------------- #    