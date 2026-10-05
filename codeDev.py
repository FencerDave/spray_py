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
    dropR = 300 * 10 ** -6  # Radius, meters (eg 300 microns)
    dRadius = 10.0**-6      # Model Solving dR Pixel Unit (1 Micron?)
    flowRate = 100.0        # g per minute total in dryer, for temperature profile,
    dropVelocity = 10.0     # m/sec initial velocity (also for thermal profile? IDK)
    
    
# REAL CLASS therefore needs _init_ and other functions?     
class SprayDroplet: # OR Spray-Arm, if not wrapping into the dryer?
    # Actually set up the droplet that we will study! 
    # HMM maybe this should be in the Spray Dryer section since these are 
    #   actually parameters of the Spraying Mechanics! 
    # THE DROPLET should be the Complicated Class that has all of the model 
    #             parameters inside of it. (So maybe I pass the dR to it? )
    def __init__(self, Slurry: SlurryData, Dryer=SprayDryer, dR:float=10**-6):
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
        self.a_lower = 4*np.pi()*self.r_vals**2 # Heat Flux LOWER
        self.a_upper = 4*np.pi()*self.r_vals**2 # Heat Flux UPPER
        self.v_voxel = (4/3*np.pi()*(self.r_vals+dR)**3) - (4/3*np.pi()*self.r_vals**3)
        #       Voxel volume is V of outer shell minus V of inner shell
        self.m_voxel = self.v_voxel * Slurry.wet_rho # g = m3 *g/m3
        #---CHEMISTRY----------------------------------------------------------
        self.solidsPct = np.ones_like(self.r_vals) * Slurry.solidsPct
        self.isDry = np.round(self.solidsPct,0)     # Binary for Solids==1
        self.k_lower = self.a_lower * Slurry.wet_condK   # W/k = W/m2K * m2(Lower)
        self.k_upper = self.a_upper * Slurry.wet_condK   # W/k = W/m2K * m2(Upper)
        self.Cp_voxel = self.m_voxel * Slurry.wet_Cp      # J/k = J/gK * g 
        self.Hvap_vox = self.m_voxel * (1-self.solidsPct) * Slurry.vapQ   
        #       J = J/g * g_WET
        #---THERMAL------------------------------------------------------------
        self.T_degC = np.ones_like(self.r_vals) * Slurry.initTemp
        self.dT_upper = np.zeros_like(self.r_vals)   # T vs r+dR (Positive)
        self.dT_lower = np.zeros_like(self.r_vals)   # T vs r-dR (Negative, 0 at center)
        self.dQdt_upper = self.dT_upper * self.k_upper    # W = W/k * dK(upper)
        self.dQdt_lower = self.dT_lower * self.k_lower    # W = W/k * dK(lower)
        self.dQdt_net = self.dQdt_upper + self.dQdt_lower # (W = W_in - W_out)
        self.dQdt_dTemp = np.zeros_like(self.r_vals)      # Eventually, residual of HVap
        self.deltaT = self.dQdt_dTemp / self.Cp_voxel     # dT = Q/Cp | dK = J / (J/K)
        return self
        
    def iterate(self,ConvectionWatts, dtime:float=10**-7):
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
        dQ_outer = dtime * ConvectionWatts # Watts * seconds = Joules
        
        # Update dT Upper and Lower values.
        for i in range(len(self.r_vals)):
            if i==0: #CASE Central Data point. dV lower = 0
                self.dT_upper[i] = self.T_degC[i] - self.T_degC[i+1]
                self.dT_lower[i] = 0 # Sphere Symmatry Boundary Cond. 
            elif i==len(self.r_vals): #CASE External Data point
                self.dT_upper[i] = 0 # Will Overwrite with Convective heat?
                self.dT_lower[i] = self.T_degC[i-1] - self.T_degC[i]
            else:
                self.dT_upper[i] = self.T_degC[i] - self.T_degC[i+1]
    
        # Calculate Heat Flux
        self.dQdt_upper = self.dT_upper * self.k_upper    # W = W/k * dK(upper)
        self.dQdt_lower = self.dT_lower * self.k_lower    # W = W/k * dK(lower)
        self.dQdt_net = self.dQdt_upper + self.dQdt_lower # (W = W_in - W_out)
        #Assign net Q to Vaporization and/or Thermal rise
        for i in range(len(self.r_vals)):
            """
            Check Solids content and Temperature. (Discuss edge-case of time step with BOTH dT and boiling?)
             * Solids >= 1 : DRY. ALL Q goes to dT
             * Temp < TVap : Wet, too cold. ALL Q goes to dT
             * Temp>=TVap and Solids <1 : BOILING. All Q goes to Evaporation. 
             * CURRENTLY NOT DEALING with edge case of T reaches TVap and excess goes to Qvap this step. 
               # IF T MEANINGFULLY OVERSHOOTS, SOLVER IS TOO LOW ACCURACY ANYWAY!!! (in fact this is a useful flag...)
            STOPPED HERE AFTER LUNCH 2026-1005
            """
        
            self.deltaT = self.dQdt_dTemp / self.Cp_voxel     # dT = Q/Cp | dK = J / (J/K)
        
        
                
        
    

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
        
    # Now Initialize/Calculate  each volume slices' enthalpy of vaporization.
    # Wonder why i don't do this for the heat capacity of the volume slice? 
    # I SHOULD, THIS IS GOOD AND THE FUTURE WAY TO HAVE THEM BE Variables, TOO! 
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
