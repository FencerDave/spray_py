# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:14:28 2026

@author: dhooke
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close("all") #Get rid of plots from last run of program
#import scipy as sp
from dataclasses import dataclass

#Use this to report better values
from math import log10, floor
def round_2(x): return round(x, 2-int(floor(log10(abs(x)))))


"""
DATA CLASSES TO INITIALIZE
 * Universe     : Constants to use? Ideal gas law? IDK
 * InletAir     : Properties, mostly temperature and humidity of the inlet airstream
 * SlurryData   : solids content, heat capacity, etc of the input slurry
 * SprayDryer   : geometry of the spray dryer, determiness temperature profile
 * SprayDroplet : initial droplet size + velocity... Anything else?
 
 Q: Where to put the dtime and dradius modeling parameters? 
     too small and it takes too long to solve, too big and models fuck up. 
 
"""
@dataclass 

class InletAir: 
    inletTemp: float = 300.0  # Temp, *C--> K
    heatTransferH: float = 1000 # Watts / m2 k Convective Heat Transfer 
    # BELOW USED TO CALCULATE THERMAL PROFILES. NOT USED YET?
    flowRate: float = 1.0     # Cubic meters per minute. NOT USED YET??
        # Medium sized spraydriers use 0.01-0.4kg/sec which is 0.5-20 m3/min
        # assume mass is 
    moisturePct: float = 0.00   # Moisture content (1.0 = 100%) - not humidity
        # Can add humidity for future physics options to help calculate.
        # This can determine heat transfer things... NOT USED YET
    heatCapacity: float = 1.005 # J / g K (Also not used yet? Will determine outlet condition...)

class SlurryData: 
    # (mostly?) NOT CURRENTLY USED: Chemical Information... Future work to use
    #                   recipe to calculate/estimate the rest of the parameters
    solvent: str = "H2O"
    solids: str = "NaCl"
    solidsPct: float = 0.00 # Solids Fraction - default is just water maybe? 
    initTemp: 25.0          # degC inlet liquid temperature
    # Physics of droplet under the "Wet" condition
    wet_condK: float = 0.5  # Watts*m/m2 K (Thermal conductivity, lookup default water)
    wet_Cp: float = 4.18    # J/g K (Heat Capacity, lookup default water)
    wet_rho: float = 1*10**6 #g/m3 (Density, default of water)
    # Physics of the REACTION PHASES
    vapTemp: float = 100.0 # degrees C for boiling the water... 
    vapQ: float = 2256.4 # j/g to evaporate water (note, only moisture% has to evap)
    # Physics of the droplet under "Dry" conditions
    dry_condK: float = 2.5 # Watts*m/m2 K... HIGHLY VARIABLE
    dry_Cp: float = 0.80    # J/g K (Heat capacity of salt or ceramic)
    dry_rho: float = 4.0*10**6 # Porous solids density... Actual is less by solidsPCT (Porosity)
    dry_damageTemp: float = 175.0 # DegC TEMPERATURE WHERE THINGS GET DAMAGED (for pass/fail criteria)
    dry_internalTemp: float = 125.0 # degC TEMPERATURE GOAL for R=0 to determine "fully dried" (for pass/fail criteria)
    

    
class SprayDryer:
    # Mostly unused for physical geometry, but this all sets up the residence
    #       time and temperature profile in the model... TBD...
    # Again start off with the currently-unused physical parameters and then 
    # get to the actual ones that determine the model for right now...
    cyl_H: float = 2.0 # m tall 
    cyl_D: float = 1.0 # m wide
    coneAngle: float = 60.0 # Degrees from Horizontal wall angle
    swirlNumber: float = 2  #No Joke they call it the "Swirl Number" - the ratio
                            # of angualar momentum to linear momentum
                            # oh this is locally calculable everywhere(?) but
                            # there's also got to be an average swirliness./
    # okay now the ones currently in use...
    dwellTime: float = 2.0  # seconds to model before the particle has left
    #                           the dryer and will no longer heat up. 
    dtime: float = 10.0**-5  # seconds per Q Solve... currently 0.01 mSecond?
    #                           (This is a model param not mechanical LOL) 
    sprayType: str = "Constant" # FLOW TYPE. THIS IMPACTS THE TEMPERATURE PROFILE
            # SprayType Options are Constant, CoCurrent, CounterCurrent?
            # CONSTANT just simplifies and keeps the Air Temp CONSTANT until we
            # develop better profiles. 
    
    
    
class SprayDroplet:
    # Actually set up the droplet that we will study! 
    dropR = 300 * 10 ** -6  # Radius, meters (eg 300 microns)
    dRadius = 10.0**-6      # Model Solving dR Pixel Unit (1 Micron?)
    flowRate = 0.1          # kg per minute total in dryer, for temperature profile,
    dropVelocity = 10.0     # m/sec initial velocity (also for thermal profile? IDK)
    
    



def sprayModel(Air, Slurry, Dryer, Drop):
    """
    Take a given droplet with its recipe and spray it through the spray Dryer
    into the inlet airstream 

    Parameters
    ----------
    inletAir : TYPE
        DESCRIPTION.
    slurryData : TYPE
        DESCRIPTION.
    sprayDryer : TYPE
        DESCRIPTION.
    sprayDroplet : TYPE
        DESCRIPTION.

    Returns
    -------
    *Series of dataframes or other matrices? where ROWS are the dR radius 
    of the drops, and the COLUMNS are the timestamps 
    
    Lifetime_Temps: Dataframe of Temperature through life
    
    Lifetime_MC:    Dataframe of Moisture content through droplet over life
    
    Results:        Dictionary? of Pass/Fail and other output recommendations 
                    for the model run. (Did we overstep TMax or Understep TMin)

    """
  
    
    # Biot Number: Geometric Ratio of Heat transfer on sphere Surface to Conduction inside.
    #      Higher Biot numbers mean slower internal conductivity, worse internal gradients.
    Biot=(Drop.dropR/3)*Air.heatTransferH/Slurry.wet_condK
    
    
    # initialize Geometry from droplet radius
    r_vals = np.arange(0, Drop.dropR + Drop.dRadius, Drop.dRadius)
    V_Vals=r_vals*0
    for i,r in enumerate(r_vals[1:]):
        V_Vals[i+1] = 4/3*np.pi*r**3 - 4/3*np.pi*r_vals[i]**3 #Volume of the shell at r - 
        
    # Now calculate each volume slices' enthalpy of vaporization.
    # Wonder why i don't do this for the heat capacity of the volume slice? 
    H_Vals = V_Vals * Slurry.VapQ * Slurry.wet_rho * (1-Slurry.solidsPct)   #Joules for each volume disc to evaporate 
    
    Temps = r_vals*0 + Slurry.initTemp #Set initial Temperature profile = T0 is constant?
    dTdt = r_vals*0 #Initialize thermal flux?  
    
    
    # NOT SURE why we're plotting things already?! 
    fig, ax = plt.subplots(figsize=[10,4])
    ax.set_ylim([0,200])
    ax.set_xlim([0,300])
    ax.set_ylabel("Temperature *C")
    ax.set_xlabel("Radius (Microns)")
    
    Plot_Time=np.arange(0,20,0.01) #Plot n solutions in 0.5 seconds in equal timesteps!
    Plots,a=np.meshgrid(r_vals, np.arange(0, Dryer.dwellTime, Dryer.dtime)) #This is how you get the data out of the loop!
    
    #------------------------------------------------------------------------
    #Okay to keep this math running, take variables back out of class for now
    Ti = Slurry.inletTemp
    k =  Slurry.wet_condK
    rho = Slurry.wet_rho
    Cp = Slurry.wet_Cp
    dt = Dryer.dtime
    dr = Drop.dRadius
    h = Air.heatTransferH
    #------------------------------------------------------------------------
    
    p=0 # Index of temperatures to plot.
    tcount=0 # index of time steps (loop until this reaches the end of time)
    t=10**-15 # Initiate Time as nonzero (avoids some divide by zero errors?) 
    
    T_S= Temps[-1] # Temperature at end of list aka @ Interface
    Ti_Variable = Ti #Initialize changing Temperature
    
    while T_S < Ti*0.95 and Temps[1]<125+273.16 and tcount<len(Plots):
        
        Ti_Variable+=(-abs((4*np.pi*(r**2))*(4*h)*(Ti_Variable-Temps[i])*dt)) #Looses by h/10? Rough Model parameter.
        Plots[tcount,:]=Temps-273.16
        tcount+=1
        #while Temps[1]<100+273.15: #Do until 
        T_S= Temps[-1] # Surface Temperature
    
        if t > Plot_Time[p]:
            p=p+1
            ax.plot(r_vals*10**6,Temps-273.16,label="time=" +str(round_2(t)*10**9)+"ns")
    
            #print(str(round_2(t)*10**3)+" milliSeconds")
    
        for i,r in enumerate(r_vals): #Loop through all r values to find dT/dt
            #Calculate Q_In, Q_Out, in units of Watts (Joules per Second.)
            #We are now trying constant time-steps again, so 
            if r==0:
                Q_out=0         #Boundary condition at middle - no place for that temperature to go.
            else:
                Q_out=(4*np.pi*(r)**2)*k*(Temps[i]-Temps[i-1])/dr    #Conduction through shell at r
    
            if r==max(r_vals):
                Q_in = (4*np.pi*(r**2))*h*(Ti_Variable-Temps[i])    #Convection at edge of shell
                
            else:
                Q_in = (4*np.pi*((r+dr)**2))*k*(Temps[i+1]-Temps[i])/dr #Conduciton through shell at -r
    
            Q = Q_in - Q_out
    
            V=4/3*np.pi*((r+dr)**3-(r)**3)
            if Temps[i]>373.16 and H_Vals[i]>0: #IF boiling point and there is Enthalpy remaining to lose:
                H_Vals[i]+= (-Q)*dt #Subtract heat from the H Vap reserve
                dTdt[i]=0 #No temp change during evaporation
            else:
                dTdt[i]= (Q) / (rho*V*Cp) # K/s,   (Watts) / (g/m3 * m3 * J/gK)
            
            
        if dTdt[-1]>0:
            t = t + dt #Capture the time elapsed
            
        #small delta T was required for smooth startup. I think we can use faster dT steps once we're in the weeds!
    #         if T_S>300:   # Once we pass 27 *C (Almost immediately, but still after ~100 loops)
    #             dt=0.005
    #         elif T_S>400: # Once we pass 127*C 
    #             dt=0.01
    
        Temps = Temps + dt * dTdt # Update the temperatures with the new temperature they have recieved?
        
    while p<len(Plots):    
        Plots[p,:]=Temps-273.15 #Grab final curve and fill it for the rest
        p+=1
    
    
    
    ax.plot(r_vals*10**6,Temps-273.15,label="time=" +str(round_2(t)*10**9)+"ns")
    ax.set_title("Heating with T="+str(round(Ti-273.16))+'*C   Biot='+str(round_2(Biot)) + '   TIME = ' + str(round_2(t))+" Seconds")
    Plots_To_Output_T.append(Plots) #Initialize list of models
    Titles_To_Output_T.append("Heating with T="+str(round(Ti-273.16))+'*C   Biot='+str(round_2(Biot)) + '   TIME = ' + str(round_2(t))+" Seconds")
    Rs_To_Output_T.append(r_vals)
    if (Temps[-1]-273.15)>175:        
        ax.axhspan(125, 175, alpha=0.1, color='red')
    else:
        ax.axhspan(125, 175, alpha=0.1, color='green')
