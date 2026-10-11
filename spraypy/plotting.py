
"""
Plotting and Visualization Functions to use in exploring Spray_Py 
and especially for use in the eventual StreamLit App! 

NOTE: SPYDER ONLY RESPECTS MATPLOTLIB (or otherwise only likes plots as png images)
 ^ And also pandas plotting or other things based on matplotlib...
 
But for Streamlit and better other data practices we should use Plotly!
SO this will contain plotting features in both libraries. 


DHooke
Initial move (separating from UnitOps.py) 2026-10-10
"""
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from spraypy.UnitOps import SprayDroplet, DropletHistory 
    # We're plotting SprayDryer info but mostly the DropletHistory result parameters. 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Recommended functions so that this code understands the type (and autocomplete)
#       for Classes made in UnitOps. 
 
"""
MOST Plotting will be done here with the DropletHistory class from SprayPy... 
Copying the definition for the class here for reference: 
    
    class DropletHistory:
        r: np.ndarray               # 1D radius axis for droplet (COLS in other arrays)
        t: np.ndarray               # 1D actual time axis (ROWS in other arrays)
        TdegC_profile: np.ndarray   # 2D radial Temperature vs time
        Solid_profile: np.ndarray   # 2D radial Solids %Pct vs time
        TdegC_avg: np.ndarray       # 1D History of Average Temperature
        Solid_avg: np.ndarray       # 1D History of 
        # ^Removed the "Surface" and "Midpoint" data here... 
        #       It's just the [0] and [-1] radius levels in the existing arrays. 
        #       A future plotter may also track the 0.25 and 0.75R points, don't make unneeded ones
        Flux_profile: np.ndarray    # 2D radial Heat flux (Joules) IN from Upper Surface
        dQdt_profile: np.ndarray    # 2D radial Net Heat change per time step
        Qnet_profile: np.ndarray    # 2D radial Net Heat Accounting (Integral of dQdt)
        # droplet history summary data points of interest
        n_records: int
        dt_min_used: int 
    
    
"""
 
def mpl_radial_plot(data:DropletHistory, t_max, t_min):
    """    
    The Default Radius-Based plot (Radius X-Axis) will have each Temperature
    records' Radius-vs-T profile, and the Solids% Profile so that we can  
    
    """
    fig,[T_ax,S_ax] = plt.subplots(2,1,sharex=True)
    fig.suptitle("Radial Profiles for T(*C) and Solids(%)")
    T_ax.set_ylabel("Temp (*C)")
    S_ax.set_ylabel("% Solids(vol)")
    S_ax.set_xlabel("Radius (um)")
    for i in range(len(data.t)):
        #Loop through all Time records
        T_ax.plot(data.r*10**6, data.TdegC_profile[i,:], 
                  label="{} sec".format(round(data.t[i],2)))
        S_ax.plot(data.r*10**6, data.Solid_profile[i,:])
    T_ax.legend()
    
    return fig
        
    
    
'''
t_array.append(theDroplet.Age)
arr_T_avg.append(theDroplet.T_degC_avg)
arr_T_outer.append(theDroplet.T_degC[-1])
arr_T_midpt.append(theDroplet.T_degC[0])
arr_Slds_avg.append(theDroplet.solidsPct_avg)
arr_Slds_outer.append(theDroplet.solidsPct[-1])
arr_Slds_midpt.append(theDroplet.solidsPct[0])

arr_H_convec.append(dQdt_convection)
arr_H_netOuter.append(theDroplet.dQdt_net[-1])
'''

#print("Is it burnt? {}".format(theDroplet.isBurnt))
#print("Is it fully dry? {}".format(theDroplet.isHappy))

#fig2 = px.line(x=theDroplet.r_vals, y=theDroplet.T_degC)
'''
fig.suptitle("Thermal Profile across Raidus")
ax[0].plot(theDroplet.r_vals*10**6, theDroplet.T_degC)
ax[1].plot(theDroplet.r_vals*10**6, theDroplet.solidsPct)
#ax[0].set_ylim(0,250)
#ax[1].set_ylim(0,1.1)
ax[0].set_ylabel("T (*C)")
ax[1].set_ylabel("Solids %")
ax[1].set_xlabel("Radius, um")

fig2, ax2 = plt.subplots(3,1,sharex=True)                 
ax2[0].plot(t_array, arr_T_avg, '-xb')                
ax2[0].plot(t_array, arr_T_outer, '--or')             
ax2[0].plot(t_array, arr_T_midpt, '-*k') 
ax2[1].plot(t_array, arr_Slds_avg, '-xb')     
ax2[1].plot(t_array, arr_Slds_outer, '--or')     
ax2[1].plot(t_array, arr_Slds_midpt, '-*k')             
ax2[2].plot(t_array, arr_H_convec, '-*y') 
ax2[2].plot(t_array,  arr_H_netOuter, '--or') 
ax2[0].set_ylabel("T (*C)")
ax2[1].set_ylabel("Avg Solids%")
ax2[2].set_ylabel("H(Watts)")
ax2[2].set_xlabel("Time(sec)")

'''

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    