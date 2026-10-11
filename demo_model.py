# -*- coding: utf-8 -*-
"""
Created on Sat 2026-10-10

Spray_Py package testing and operations, WIP to remove any module-level scripts
so that it can act as a package and be imported by pytest and also to be built
into streamlit apps and other user-focused sides of the proejct. 

Model is working enough when run manually with the self.iterate() command and
given a heat value per each timestep... Enough to reach this stopping point 
and functionalize / deploy this as a little package. 

@author: dhooke
"""
from spraypy import InletAir, SlurryData, SprayDryer, SprayDroplet
# Effectively, import all... 
from spraypy.plotting import mpl_radial_plot as plt_rad


theAir = InletAir()
theSlurry = SlurryData()                
theDryer = SprayDryer(dropR= 125*10**-6)

#defaultDrop = SprayDroplet(theSlurry, theDryer, dR=10**-6)
#defaultDrop.iterate(0)
theDroplet = SprayDroplet(theSlurry, theDryer, dR=1*10**-6)
                
#import plotly.express as px
import matplotlib.pyplot as plt
import time 
#import plotly.graph_objects as go

dt = 100*10**-9             # Model Delta-Time in sec (Nanosec so long as *10**-9)
#                           # For reference, Our Dyamic dt minimum limits 
                            #     WET @ h=1000:  2000 *10**-9?
                            #     DRY Worst-Case? : 53 * 10**-9?
#N_Datapoints = 50*10**6     # millions of points. 
t_drying = 0.25             # time in seconds
t_DYNAMIC_Q = True 


N_Datapoints = round(t_drying/dt)
N_Readouts = 20
t_Trigger_step = t_drying / N_Readouts # milisecs to data log
t_Trigger_next = 0
t_max = dt*N_Datapoints 

t_array = []
arr_T_avg   = []
arr_T_outer = []
arr_T_midpt   = []

arr_Slds_avg    = []
arr_Slds_outer    = []
arr_Slds_midpt    = []

arr_H_convec = []
arr_H_netOuter = []

# Loop through 2 seconds and every 0.1 seconds update the plot
# Watts = k * W/m2K * m2
i = 0
fig, ax = plt.subplots(2,1, sharex=True)
#ax[0].plot(theDroplet.r_vals*10**6, theDroplet.T_degC)
#ax[1].plot(theDroplet.r_vals*10**6, theDroplet.solidsPct)


time_100 = 0 #Initialize clock readout...
time_est = -666*60 #initialize... Nonsense Readings are a symbol
while theDroplet.Age < t_drying: # Fix to track Dynamic "TRUE" Timing! 
    
    # As soon as we have at iterated throuh N= Number of voxels times,
    # Time the function and use that to predict our model run time
    R_points = len(theDroplet.r_vals)
    if i == round(R_points):
        time_100 = time.perf_counter()
        # RESET Timer to readout... 
    elif i == round (R_points + 100):
        time_per =  (time.perf_counter() - time_100)/100.0 #Avg 100 loops?
        time_est = N_Datapoints * time_per
        print ("Running {} loops... Should take {} min...".format(
            N_Datapoints, round((time_est)/60)))
        # NOTE: Still predicting RunTime from the time-steps used at the BEGINNING
    else:
        pass
    
    dTemp_outer = theAir.inletTemp - theDroplet.T_degC[-1]
    dQdt_convection = dTemp_outer * theAir.heatTransferH * theDroplet.a_upper[-1]
    
    #theDroplet.iterate(dQdt_convection, dt, dynamic_dt = True)
    theDroplet.iterate(dQdt_convection, dt, dynamic_dt = t_DYNAMIC_Q)
    
    i=i+1
    #if round(theDroplet.Age*10**9) % (t_Trigger_msec*10**6) == 0 : 
    if theDroplet.Age > t_Trigger_next: 
        t_Trigger_next += t_Trigger_step
            
        ax[0].plot(theDroplet.r_vals*10**6, theDroplet.T_degC)
        ax[1].plot(theDroplet.r_vals*10**6, theDroplet.solidsPct)
        #ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_lower)
        #ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_upper)
        print("time: {}/{} m-sec. TMax = {}  EnergyCheck: {}%  Elapsed: {}/{} min".format(
            round(theDroplet.Age*10**3),round(t_max*10**3), round(max(theDroplet.T_degC)), 
            round(100*sum(theDroplet.Qnet_vox)/theDroplet.Qnet_drop),
            round((time.perf_counter() - time_100)/60), round(time_est/60)))
        
        t_array.append(theDroplet.Age)
        arr_T_avg.append(theDroplet.T_degC_avg)
        arr_T_outer.append(theDroplet.T_degC[-1])
        arr_T_midpt.append(theDroplet.T_degC[0])
        arr_Slds_avg.append(theDroplet.solidsPct_avg)
        arr_Slds_outer.append(theDroplet.solidsPct[-1])
        arr_Slds_midpt.append(theDroplet.solidsPct[0])
        
        arr_H_convec.append(dQdt_convection)
        arr_H_netOuter.append(theDroplet.dQdt_net[-1])
        

#print("Is it burnt? {}".format(theDroplet.isBurnt))
#print("Is it fully dry? {}".format(theDroplet.isHappy))

#fig2 = px.line(x=theDroplet.r_vals, y=theDroplet.T_degC)
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