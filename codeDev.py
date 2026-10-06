# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:14:28 2026

@author: dhooke
"""
# %reset -f

import numpy as np
import matplotlib.pyplot as plt
plt.close("all") #Get rid of plots from last run of program
#import scipy as sp
from dataclasses import dataclass

#Use this to report better values
from math import log10, floor
def round_2(x): return round(x, 2-int(floor(log10(abs(x)))))


"""
NOTE ON UNITS: 
    USING M-G-S Units, so METERS for length scales (not sure why but that's what I have)
                          GRAMS for mass scales
                          JOULES for heat
                          SECONDS for time
                          CELCIUS for temperature (checked. No need for Kelvin)
                          


DATA CLASSES TO INITIALIZE
 * Universe     : Constants to use? Ideal gas law? IDK
 * InletAir     : Properties, mostly temperature and humidity of the inlet airstream
 * SlurryData   : solids content, heat capacity, etc of the input slurry
 * SprayDryer   : geometry of the spray dryer, determiness temperature profile
                     # ALSO initial droplet size + velocity... The Spray Arm
NORMAL CLASS: SprayDroplet
 * __init__     : Initialize droplet... With or without pre-determined history?
 * iterate      : Take an external heat(power) value and dtime step and update
                     # across the entire droplet one step...
 * setup_history: Optional set up to initialize the history tables? 
 * record_history:Record the current droplet state in the history dataframe(s)?
 
 
 
 Q: Where to put the dtime and dradius modeling parameters? 
     too small and it takes too long to solve, too big and models fuck up. 
 
"""
@dataclass 
class InletAir: 
    inletTemp: float = 300.0  # Temp, *C--> K
    heatTransferH: float = 1000 # Watts / m2 k Convective Heat Transfer
                    # I think it's normally 1000 but changing to get it to work
    # BELOW USED TO CALCULATE THERMAL PROFILES. NOT USED YET?
    flowRate: float = 1.0     # Cubic meters per minute. NOT USED YET??
        # Medium sized spraydriers use 0.01-0.4kg/sec which is 0.5-20 m3/min
        # assume mass is 
    moisturePct: float = 0.00   # Moisture content (1.0 = 100%) - not humidity
        # Can add humidity for future physics options to help calculate.
        # This can determine heat transfer things... NOT USED YET
    heatCapacity: float = 1.005 # J / g K (Also not used yet? Will determine outlet condition...)

@dataclass 
class SlurryData: 
    # (mostly?) NOT CURRENTLY USED: Chemical Information... Future work to use
    #                   recipe to calculate/estimate the rest of the parameters
    solvent: str = "H2O"
    solids: str = "NaCl"
    solidsPct: float = 0.5 # Solids Fraction - default is just water maybe? 
    initTemp: float = 25.0          # degC inlet liquid temperature
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
    dry_damageTemp: float = 175.0 # DegC TEMPERATURE WHERE THINGS GET DAMAGED (for pass/fail criteria)
    dry_targetTemp: float = 125.0 # degC TEMPERATURE GOAL for R=0 to determine "fully dried" (for pass/fail criteria)
    

@dataclass 
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
    #                           OPTIONAL NAN? Keep going until droplet reaches 
    #                               internal T = 95% of External Temp? (Until
    #                               other definition of "steady state"?) 
    dtime: float = 10.0**-5  # seconds per Q Solve... currently 0.01 mSecond?
    #                           (This is a model param not mechanical LOL) 
    sprayType: str = "Constant" # FLOW TYPE. THIS IMPACTS THE TEMPERATURE PROFILE
            # SprayType Options are Constant, CoCurrent, CounterCurrent?
            # CONSTANT just simplifies and keeps the Air Temp CONSTANT until we
            # develop better profiles. 
    dropR: float = 300 * 10 ** -6  # Radius, meters (eg 300 microns)
    dRadius: float = 10.0**-6      # Model Solving dR Pixel Unit (1 Micron?)
    flowRate: float = 100.0        # g per minute total in dryer, for temperature profile,
    dropVelocity: float = 10.0     # m/sec initial velocity (also for thermal profile? IDK)
    
    
# REAL CLASS therefore needs _init_ and other functions?     
class SprayDroplet: 
    '''
    SprayDroplet: The Model itself! 
    
    '''
    def __init__(self, Slurry: SlurryData, Dryer: SprayDryer, dR:float=1*10**-6):
        """
        Initialize a Droplet
         * INPUT slurry parameters (wet physical properties, solidsPct, initial Temp)
         * INPUT droplet size (Eventually from Sprayer Setup)
         * INPUT droplet dRadius for model solve 
         #          (Default 0.1um OR default to 0.1% of radius for 1k pts??)
         
         * Output X-Axis (Raidius Array)
         * OUTPUT initial Temperature Array
         * Output initial Moisture Content Array (or solids content damnit... Decide)
         * Output Zeroes dT slope Array (future Q: Avg or Right-side?)
                                     # EULER or Reverse-Euler or TriDiagonal??
         * Output Wet (initial) Inlet/Outlet Conductivity arrays? (of each Voxel, J/k)
                                     # Either explicit or as one (In/Out differ by dr)
         * Output Wet (initial) Heat Capacity array (of each Voxel, J/k)
         * Output Wet (initial) Heat of Vapor array (of each Voxel, J per _%mc?_ Unit?)
         * Output WARNINGS / TARGETS (inherited from Slurry Recipe?)
         OUTPUT AS A CLASS OF DROPLET with THOSE ARRAYS... (Dataframe or just list?)
        """
        #---GEOMETRY-----------------------------------------------------------
        self.r_vals = np.arange(0, Dryer.dropR - dR, dR) # Radius Axis
        #   Each VOXEL is idenfified above then by its LOWER r value.
        #       Thus the voxel area and volume calcs are from R to R+dR
        #       This is why the maximum r is (R-dR)... REMEMBER THAT! 
        self.a_lower = 4*np.pi*self.r_vals**2 # Heat Flux LOWER
        self.a_upper = 4*np.pi*(self.r_vals + dR)**2 # Heat Flux UPPER
        self.v_voxel = (4/3*np.pi*(self.r_vals+dR)**3) - (4/3*np.pi*self.r_vals**3)
        #       Voxel volume is V of outer shell minus V of inner shell
        self.m_voxel = self.v_voxel * Slurry.wet_rho # g = m3 *g/m3
        self.m_solid = self.m_voxel * Slurry.solidsPct 
        self.m_liquid = self.m_voxel * (1-Slurry.solidsPct)
        #---CHEMISTRY----------------------------------------------------------
        self.solidsPct = np.ones_like(self.r_vals) * Slurry.solidsPct
        self.isDry = np.floor(self.solidsPct)     # Binary for Solids==1
        self.k_lower = self.a_lower / dR * Slurry.wet_condK   # W/k = Wm/m2K * m2(Lower) /m_dR
        self.k_upper = self.a_upper / dR * Slurry.wet_condK   # W/k = Wm/m2K * m2(Upper) /m_dR
        self.Cp_voxel_wet = self.m_voxel * Slurry.wet_Cp      # J/k = J/gK * g 
        self.Cp_voxel_dry = self.m_solid * Slurry.dry_Cp      # J/k = J/gK * g 
        self.Hvap_vox = self.m_liquid * Slurry.vapQ   # J = g_liquid * J/gLiqVap
        self.H_dT_vox = self.Cp_voxel_wet * (Slurry.vapTemp - Slurry.initTemp) # J = J/k * dK
        # Copying Slurry Constants for future use, including all dry params...
        self.solidsPct_init = Slurry.solidsPct
        self.initTemp       = Slurry.initTemp
        self.vapQ           = Slurry.vapQ
        self.vapTemp        = Slurry.vapTemp
        self.dry_condK   = Slurry.dry_condK    # (remember to re-divide by dR because its Wm/m2K)
        self.wet_condK   = Slurry.wet_condK
        self.dry_damageTemp = Slurry.dry_damageTemp
        self.dry_targetTemp = Slurry.dry_targetTemp
        #       J = J/g * g_WET
        #---THERMAL------------------------------------------------------------
        self.T_degC = np.ones_like(self.r_vals) * Slurry.initTemp
        self.dT_upper = np.zeros_like(self.r_vals)   # T vs r+dR (Positive)
        self.dT_lower = np.zeros_like(self.r_vals)   # T vs r-dR (Negative, 0 at center)
        self.dQdt_upper = self.dT_upper * self.k_upper    # W = W/k * dK(upper)
        self.dQdt_lower = self.dT_lower * self.k_lower    # W = W/k * dK(lower)
        self.dQdt_net = self.dQdt_upper + self.dQdt_lower # (W = W_in - W_out)
        self.Qnet_vox = np.zeros_like(self.r_vals) # TRACK THE TOTAL HEAT GAINED OVER TIME
        # STRATEGY post Claude-consult: BOOK-KEEP focuses on HEAT and then updates
        #           Temperature & Moisture based on the current Energy in the Voxel
        #           so now, DON'T TOUCH Hvap and H_dT, but COMPARE those to Qnet for the voxel
        self.dR = dR       
        
    def iterate(self,ConvectionWatts, dtime:float=10**-7, debug=False):
        """
        Iterate a Droplet:
         * INPUT Droplet Parameters (Listed Above)
         * INPUT dtime (turns Watts into Joules)
         * INPUT External Heat (Watt or Joules) from Convection?
         * CALCULATE Inlet Heat (Cond or External per voxel)
         * CALCULATE Outlet Heat (Cond or Zero per voxel)
         * CALCUALTE Moisture Loss (If TBoil is reached and not Dry)
         * CALCULATE Temperature Rise (If not TBoil or YesDry)
                 # EVENTUALLY move away from Forward-EULER? 
                 # EVENTUAL ACCURACY - Evap and dT could be smoothed towards reality
         * UPDATE Temperature, Moisture Content
         * UPDATE Conductivity, Heat Capacity based on new Moisture 
         * OUTPUT Droplet Parameters
         * OUTPUT Warnings/ Target Flags (Over TMax, or Reached T-target@R0, etc?)
        """
                
        # Update dT Upper and Lower values.
        # Better or worse than the literal gradient to each side?? hmm... 
        # MAY lead to some Instsability... Come back to this. 
        self.dT_upper = np.gradient(self.T_degC)
        self.dT_lower = -np.gradient(self.T_degC)
        '''
        for i in range(len(self.r_vals)):
            # Positive Slope for Upper (Gaining Heat)
            # Negative slope for Lower (Losing Heat)
            # IDEA: Replace ALL OF THIS with np.gradient()
            #           (One gradient per cell, dQdt difference depends on area In/out difference)
            if i==0: #CASE Central Data point. dV lower = 0
                self.dT_upper[i] = self.T_degC[i+1] - self.T_degC[i]
                self.dT_lower[i] = 0 # Sphere Symmatry Boundary Cond. 
                
            elif i==len(self.r_vals)-1: #CASE External Data point
                self.dT_upper[i] = 0 # Will Overwrite with Convective heat?
                self.dT_lower[i] = self.T_degC[i-1] - self.T_degC[i]
                
            else:
                self.dT_upper[i] = self.T_degC[i+1] - self.T_degC[i]
                self.dT_lower[i] = self.T_degC[i-1] - self.T_degC[i]
        '''
        # Calculate Heat Flux
        self.dQdt_upper = self.dT_upper * self.k_upper    # W = W/k * dK(upper)
        self.dQdt_upper[-1] = ConvectionWatts             # Overwrite Convection dQdt
        self.dQdt_lower = self.dT_lower * self.k_lower    # W = W/k * dK(lower)
        self.dQdt_net = self.dQdt_upper + self.dQdt_lower # (W = W_in - W_out)
        
        self.Qnet_vox = self.Qnet_vox + (self.dQdt_net * dt) # Q += Watts*secs
        
        
        # BASED ON Q Net, Update TEMPERATURE!
        # Interpolate Temperature where wet based on Cp wet
        self.T_Wet = (self.Qnet_vox <= self.H_dT_vox) * (self.Qnet_vox - self.H_dT_vox) / (self.initTemp - self.vapTemp)
        self.T_Boil = (self.Qnet_vox > self.H_dT_vox) * (self.Qnet_vox <= (self.H_dT_vox + self.Hvap_vox) ) * self.vapTemp
        self.T_Dry = (self.Qnet_vox > self.Hvap_vox ) * (self.Qnet_vox-self.Hvap_vox -self.H_dT_vox)/self.Cp_voxel_dry
        self.T_degC = self.T_Wet + self.T_Boil + self.T_Dry
        
        # IF this logic works, use same boolean math (now on Temperature?) to Reassign conductivity K values?
        self.k_lower = ((self.T_degC > self.vapTemp) * (self.a_lower/self.dR * self.dry_condK)
                       +(self.T_degC <= self.vapTemp) * (self.a_lower/self.dR * self.wet_condK))
        self.k_upper = ((self.T_degC > self.vapTemp) * (self.a_upper/self.dR * self.dry_condK)
                       +(self.T_degC <= self.vapTemp) * (self.a_upper/self.dR * self.wet_condK))
        
        # Calculate Report Moisture Content based on Q Net
        self.Solids_Wet = (self.Qnet_vox <= self.H_dT_vox) * self.solidsPct_init
        self.Solids_Boil = (self.Qnet_vox > self.H_dT_vox) * (self.Qnet_vox <=(self.H_dT_vox + self.Hvap_vox))*(
                            (self.Qnet_vox - self.H_dT_vox-self.Hvap_vox)/(self.solidsPct_init - 1.00))
        self.Solids_Dry = (self.Qnet_vox > self.Hvap_vox ) * 1.0
        
        
        '''#Assign net Q to Vaporization and/or Thermal rise
        for i in range(len(self.r_vals)):
            """
            Check Solids content and Temperature. (Discuss edge-case of time step with BOTH dT and boiling?)
             * Solids >= 1 : DRY. ALL Q goes to dT
             * Temp < TVap : Wet, too cold. ALL Q goes to dT
             * Temp>=TVap and Solids <1 : BOILING. All Q goes to Evaporation. 
             * CURRENTLY NOT DEALING with edge case of T reaches TVap and excess goes to Qvap this step. 
               # IF T MEANINGFULLY OVERSHOOTS, SOLVER IS TOO LOW ACCURACY ANYWAY!!! (in fact this is a useful flag...)
               # This would be the only reason to separate "dQdt_net" into dT and dRXN parts... Ignore for now!
            
            
            *RATHER THAN LOOPING THROUGH... 
              * Make a Boolean Array for "Rows where this is true"
              * Adjust all of those rows by the amount needed? Don't Adjust anything to the others?'
            """
            if not(self.isDry[i]) and self.T_degC[i] >= self.vapTemp:              
                # IF Not Dry and Have reached Temperature: 
                # Two Scenarios: 
                    # This Heat Packet IS ENOUGH to dry it out 
                    #       (therefore it will be DRY, AND will have Q left to gain dT)
                    # This heat Packet IS NOT ENOUGH to dry it out
                    #       (therefore remove Q from the HVap counter, Recalculate Solids%)
                dQ = self.dQdt_net[i] * dtime
                if dQ >= self.Hvap_vox[i]:
                    # DRY THIS VOXEL! 
                    # Subtract this Q from NET (which will still go to heating)
                    self.isDry[i] = 1 # Affirm Now it's dry
                    self.dQdt_dTemp[i] = self.dQdt_net[i] - self.Hvap_vox[i]/dt # W = W - dJ/dsec
                    self.Hvap_vox[i] = 0        # No more liquid to evaporate
                    self.m_liquid[i] = 0        # No more liquid to evaporate
                    self.m_voxel[i] = self.m_solid[i] # No more liquid
                    # Now Re-Assign Properties to the "DRY" Chemical Cp, K
                    self.Cp_voxel[i] = self.m_voxel[i] * self.dry_Cp     # W/k = W/g * g
                    self.k_lower[i] = self.a_lower[i] * self.dry_condK_dR   # W/k = W/m2K * m2(Lower)
                    self.k_upper[i] = self.a_upper[i] * self.dry_condK_dR   # W/k = W/m2K * m2(Upper)
                else:
                    # Remove dQ amount of Hvap (and liquid) from the Voxel
                    # NO MORE dQdt_net to go to heating! 
                    # Subtract heat from Hvap and recalculate moisture%v etc
                    self.isDry[i] = 0 # Affirm not yet dry
                    self.Hvap_vox[i] = self.Hvap_vox[i] - dQ  
                    self.m_liquid[i] = self.Hvap_vox[i] / self.vapQ
                    self.m_voxel[i] = self.m_solid[i] + self.m_liquid[i]
                    self.dQdt_dTemp[i] = 0 # No more Net Heat for dTemp
                # Now adjust solidsPCT (to either 1 if dry or )
                self.solidsPct[i] = self.m_solid[i] / self.m_voxel[i]
                # OPTIONAL Here: Snap T to boiling point ()
                self.T_degC[i] = self.vapTemp # HARD SNAP to T_VAP? (Temporary fix?)
                # DONE the Boiling Case, and dQ Net now all goes to Heating
                if debug: breakpoint()
            else:
                self.dQdt_dTemp = self.dQdt_net #IF no evaporation, all Q goes to dTemp
            # END FOR-Loop looping through all r values 
           
        # NOW apply any remaining dQdt Net (Power) to become Heating
        self.deltaT = self.dQdt_dTemp * dtime / self.Cp_voxel     # dT = Q/Cp | dK = W*sec / (J/K)    
        #Now incriment any X values where Temperature has increased
        self.T_degC = np.add(self.T_degC, self.deltaT).tolist()
        #Now flag the system, for either reaching the Too-Hot external
        #       or for reaching the target internal Temperature
        '''
        self.isBurnt = bool(max(self.T_degC) >= self.dry_damageTemp)
        self.isHappy = bool(min(self.T_degC) >= self.dry_targetTemp)
            
        
# Generate Defaults
theAir = InletAir()
theSlurry = SlurryData()                
theDryer = SprayDryer()

defaultDrop = SprayDroplet(theSlurry, theDryer, dR=10**-6)

theDroplet = SprayDroplet(theSlurry, theDryer, dR=10**-6)
                
#import plotly.express as px
import matplotlib.pyplot as plt
#import timeit 
#import plotly.graph_objects as go

#fig = px.line(x=theDroplet.r_vals, y=theDroplet.T_degC)    

dt = 5*10**-7
N_Datapoints = 2*10**3
N_Readouts = 50
Trigger=round(N_Datapoints/N_Readouts)
t_max = dt*N_Datapoints
time=0
debugLoopX=3
# Loop through 2 seconds and every 0.1 seconds update the plot
# Watts = k * W/m2K * m2
i = 0
fig, ax = plt.subplots(2,1, sharex=True)
fig2, ax2 = plt.subplots(1,1, sharex=True)
ax[0].plot(theDroplet.r_vals*10**6, theDroplet.T_degC)
ax[1].plot(theDroplet.r_vals*10**6, theDroplet.solidsPct)
ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_lower)
ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_upper)

while i<N_Datapoints:
    time = time + dt
    dTemp_outer = theAir.inletTemp - theDroplet.T_degC[-1]
    dQdt_convection = dTemp_outer * theAir.heatTransferH * theDroplet.a_upper[-1]
    theDroplet.iterate(dQdt_convection, dt)
    
    i=i+1
    if i/Trigger == np.floor(i/Trigger): # True for 1000 but not for 1001 
            
        ax[0].plot(theDroplet.r_vals*10**6, theDroplet.T_degC)
        ax[1].plot(theDroplet.r_vals*10**6, theDroplet.solidsPct)
        ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_lower)
        ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_upper)
        print("time: {}/{} microsec".format(round(time*10**6),round(t_max*10**6) ))

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


ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_lower)
ax2.plot(theDroplet.r_vals*10**6, theDroplet.dT_upper)

                
#fig.show(renderer="browser")
#fig2.show(renderer="browser")
#fig3.show()
        
    

"""
LIST OF FUNCTIONS TO BUILD

Initialize a Droplet
 * INPUT slurry parameters (wet physical properties, solidsPct, initial Temp)
 * INPUT droplet size (Eventually from Sprayer Setup)
 * INPUT droplet dRadius for model solve (OR default to 0.1% of radius for 1k pts)
 * Output X-Axis (Raidius Array)
 * OUTPUT initial Temperature Array
 * Output initial Moisture Content Array (or solids content damnit... Decide)
 * Output Zeroes dT slope Array (future Q: Avg or Right-side?)
                             # EULER or Reverse-Euler or TriDiagonal??
 * Output Wet (initial) Inlet/Outlet Conductivity arrays? (of each Voxel, J/k)
                             # Either explicit or as one (In/Out differ by dr)
 * Output Wet (initial) Heat Capacity array (of each Voxel, J/k)
 * Output Wet (initial) Heat of Vapor array (of each Voxel, J per _%mc?_ Unit?)
 * Output WARNINGS / TARGETS (inherited from Slurry Recipe?)
 OUTPUT AS A CLASS OF DROPLET with THOSE ARRAYS... (Dataframe or just list?)



Determine Droplet Lifetime Thermal Profile? 
# TO INITIALIZE, CAN ALSO BE MADE VERY SIMPLIFIED: 
    1) CONSTANT TEMPERATURE (Giant Sprayer and/or well-mixed steady state)
    2) LINEAR TEMPERATURE INLET-->OUTLET (Rough Approx of Cocurrent)
    3) DUAL-LINEAR or Quadratic Inlet-->Outlet (Rough Approx of CounterSpray)
 * INPUT Overriding Assumptions (SEE ABOVE)
 * INPUT Air Properties (Temperature, Mass/Volume Flow and Heat Capacity)
 * INPUT Bulk Fluid Properties (Tinlet, Heat capaicty & Latent heat)
 * INPUT Spray Dryer Properties (Direction, Flow Rate, Particle Velocity)
 * DETERMINE if Simplified Assumptions to be used
 * CALCULATE Residence Time Estimate (Particle Drag? Velocities and Airflow???)
 * CALCULATE Thermal Balance (Net Adiabadic, What is T_Outlet?)
 * OUTPUT Lookup Table for Droplet T_External (for heat flux) vs Time))
 * OUTPUT Lookup Table for Heat Transfer Coefficients vs time?? (IDK Physics)


Determine Droplet External Heat Transfer from Convection:
 * INPUT Droplet External Temp at current Time (T @ R=RMax) 
 * INPUT Airflow Temperature at current Time (from thermal Profile)
 * INPUT Airflow Convective Coefficient at current Time? (or constant?)
 * OUTPUT Droplet Boundary Watts per m2 (NOT JOULES so that only ITERATE has dtime Component)


Iterate a Droplet:
 * INPUT Droplet Parameters (Listed Above)
 * INPUT dtime (turns Watts into Joules)
 * INPUT External Heat (Watt or Joules) from Convection?
 * CALCULATE Inlet Heat (Cond or External per voxel)
 * CALCULATE Outlet Heat (Cond or Zero per voxel)
 * CALCUALTE Moisture Loss (If TBoil is reached and not Dry)
 * CALCULATE Temperature Rise (If not TBoil or YesDry)
         # EVENTUALLY move away from Forward-EULER? 
         # EVENTUAL ACCURACY - Evap and dT could be smoothed towards reality
 * UPDATE Temperature, Moisture Content
 * UPDATE Conductivity, Heat Capacity based on new Moisture 
 * OUTPUT Droplet Parameters
 * OUTPUT Warnings/ Target Flags (Over TMax, or Reached T-target@R0, etc?)
 
Report Droplet Values of Note: 
        # What are we interested in plotting??
 * INPUT Droplet Parameters
 * OUTPUT Droplet Mass-Averaged (or Volume-Averaged?) Temperature
 * OUTPUT Droplet Mass-Averaged Moisture %? 
 * OUTPUT Surface and Internal (r=R and r=0) Temperatures
 * OUTPUT %mass Under Tmin target, and %mass Over TMax Flag??


Initialize Data Table?:
    # See notes on data logging below? 
 * Input Initialized Droplet (Therefore correct length of Radius table)
 * Input planned record Triggers/Timesteps?
 * Output 2x2 Matrices for each parameter to record?
 * For Record Keeping, INCLUDE a 2x2 Matrix (?) with TIME as value?. 

Log Data into Table:
     # A better simulation would turn this into matrix math and do this 
     #      all at once. I kind of forget how to do that so my instinct is to 
     #      just create the chart/matrix as we go? Can initialize the 2D matrices
     #      at the beginning of a model run, or else 
 * Input Droplet Parameters
 * Input Current Time
 * Input Planned Record Triggers/ Timesteps? 
 * Input Record Matrices? Initialized elsewhere.
 * Output Record Matrices with new line filled in...
 
Log Droplet Values-of-Note into Table?
 * INPUT finished Record Matrices
 * RUN "Report Droplet Values of note" on each timestep in table
 * OUTPUT Table / Dataframe with columns for Time, Tavg, Tmin, Tmax, %mc, etc.
 
 Plot Thermal History: 
 * Input Record Matrix of Temperature and Time
 * Input Target and OverTemp flags/limits
 * OUTPUT Plot of T-vs-R with lines for each time
 #          Include Legend (maybe input decisions here) for time and subTime 
 #          (such as, Bold every 10 Time-lines to count CentiSeconds?)
 * OUTPUT Tmin and TMax box on the chart
 #          COLOR by whether (A) T>TMin @ R=0 and/or (B) T>TMax @ at r=R
 #              IF A & !B, SUCCESS
 #              IF A & B, Dried but OVERHEATED - Lower Temperature
 #              IF !A & !B, UNDERDRIED but didn't Overheat - Raise Temperature
 #              IF !A & B, FAILURE - PSD and Dryer Mismatch! 
 
 Plot Values of Note vs Time
 * Input Values-of-note vs time
 * OUTPUT Temperatures of note vs time
 * OUTPUT Moisture% vs time
 * COLORIZE "Pass" and "Fail" Times if they exist, as well as box from before.
 
 
 

"""










