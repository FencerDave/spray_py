# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:14:28 2026

@author: dhooke
"""
# %reset -f

import numpy as np
#import scipy as sp

from dataclasses import dataclass
import warnings 

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
                    # Ive historically used 1000 = a super turbulent spray dryer
                    # Nussult Number"Nu=2" and other notes on "Stationary" 
                    #       spheres give 100-200 for 600um Diameter, but we are
                    #       VERY NOT STATIONARY! In any case, future physics!
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
    dry_condK: float = 2.5 # Watts*m/m2 K... HIGHLY VARIABLE per material...
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
    #FUTURE USE: Droplet Mechanics and STATISTICS 
    dropR: float = 300 * 10 ** -6  # Radius, meters (eg 300 microns)
    dropR_sigmaDist = 2 # FUTURE USE - SPREAD of droplet Size Distribution.
     # ^ Eventaully, Model a "Larger" and "Smaller" particle from distribution to 
     #               check whether the largest are un-dried or the smallest are burnt!    
     
    flowRate: float = 100.0        # g per minute total in dryer, for temperature profile,
    dropVelocity: float = 10.0     # m/sec initial velocity (also for thermal profile? IDK)
    
@dataclass
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
    # def to_pd_df(self) -> pd.DataFrame: # Not sure what this entails yet?
    
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
        self.r_vals = np.arange(0, Dryer.dropR, dR) # Radius Axis
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
        self.m_liquid_init = self.m_voxel * (1-Slurry.solidsPct)
        
        #---CHEMISTRY----------------------------------------------------------
        self.solidsPct = np.ones_like(self.r_vals) * Slurry.solidsPct
        self.isDry = np.floor(self.solidsPct)     # Binary for Solids==1
        #self.k_lower = self.a_lower / dR * Slurry.wet_condK   # W/k = Wm/m2K * m2(Lower) /m_dR
        # Eliminated k_lower in order to only have one conduction per Area, apply to above/below as needed
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
        #self.dQdt_upper = self.dT_upper * self.k_upper    # W = W/k * dK(upper)
        #self.dQdt_lower = self.dT_lower * self.k_lower    # W = W/k * dK(lower)
        #self.dQdt_net = self.dQdt_upper + self.dQdt_lower # (W = W_in - W_out)
        self.Qnet_vox = np.zeros_like(self.r_vals) # TRACK THE TOTAL HEAT GAINED OVER TIME
        self.Qnet_drop = 0.0 #Track and compare to the sum of Qnet to identify losses
        
        #---MODELING AND HISTORY-----------------------------------------------
        self.Age = 0.0       #ALSO Track the age rather than doing that outside of loops. 
                                #^Especially important with have DYNAMIC dTime 
        self.dR = dR       
        self.history: DropletHistory # Empty arrays? Not sure how to initialize? 
        
    def _stable_dtime(self, dt_check=0, margin: float=0.0, geoAlpha: float=3.0,
                      dt_min = 10**-9, dt_max = 10**-6):
        """
        Given a droplets current status, Use Biot and Fourier Math to estimate
        the current minimum delta time step for the heating model to use. 
        CALCULATE:
            Bi_drop:    float, Biot Number with L = R/3 for sphere
            Bi_grid:    float, Biot Number with L = dR for gridmesh
            
            C_sigG :    float, dt limit by VOXEL Heat capacity divided by all conductances
                        (so [J/K]/[W/K]=[sec] even though the dK terms are different)
                        (This dt will be worst at CENTRAL VOXEL always? since
                         no inner conduction AND smallest heat capacity in J/K)
                        
        
        INPUT:
            dt_check:   float OR BOOL?: OPTIONAL: IF GIVEN, check this dt value
                            against the caluclated dtime recommendation. 
                            Flag a WARNING message if dt is NOT COMPLIANT.
                            (Default == 0 will always pass and therefore be ignored)
            margin:     float, Margin by which to go beyond the theoretical 
                        literature-recommended dtime step... 
                        dTime gets DIVIDED BY (1 + MARGIN) SO FOR INSTANCE:
                        m = 0.00    Returns THEORETICAL LIMIT dTime
                        m = 1       Returns dTime(th)/2
                        m = 0.5     Returns dTime(th)/1.5
                        m = -0.1    Returns dTime(th)/0.9 # FAILS LIMIT ON PURPOSE?
        RETURN: 
            dt_step:    float, recommended dTime step based on fourier limit
                        and any additional margins given....
                        dt_lim  = dR^2 * rho * Cp / k 
                                = [m2][g/m3][J/gK] / [Wm/m2K]
                        dt_step = dt_lim / (geometry * (1 + margin))
        """
        # Calculate per-Voxel Heat Capacity (BoolMath from wet vs dry via Hnet)
        #   This way we dont use Dry values at the center until it reaches the center. 
        #   Dealing with BOILING HERE is a bit tricky? CP is infinite-ish since T won't change but also there is the Q associated with fully evaporating...? 

        self.Cp_voxel = np.where(self.Qnet_vox <= self.H_dT_vox, self.Cp_voxel_wet,
                        np.where(self.Qnet_vox < (self.H_dT_vox + self.Hvap_vox), np.inf, self.Cp_voxel_dry))
        
        sigG_vox = np.zeros_like(self.Cp_voxel)
        sigG_vox[:-1] += self.k_upper[:-1]            # interface to the outer neighbour
        sigG_vox[1:] += self.k_upper[:-1]             # interface to the inner neighbour
        sigG_vox[-1] += self.k_upper[-1] * self.a_upper[-1]    # MISSING CONVECTION TERM, but Instead double counting Conduction at edge?
            
        dt_C_SigG = self.Cp_voxel / (geoAlpha * sigG_vox) # Ratio J/W locally and divided by Geometric 01Alpha" from literature? 
        # CITE ANY SOURCE FOR THAT? -->>>>>
        
        dt_RESULT = min(dt_C_SigG) / (1 + margin) 
        self.dt_array = dt_C_SigG / (1 + margin) # TEMPORARY OUTPUT TO CHECK dtime Values??
        
        if dt_check > min(dt_C_SigG):
            warnings.warn( "dt is LARGER than Stability Limit per Euler/Fourier",
              category=RuntimeWarning, stacklevel=2)
        elif dt_check >= (min(dt_C_SigG)/(1+margin)):
            warnings.warn( "dt is Within MARGIN of theoretical Limit Eu/Bi/Fo",
              category=RuntimeWarning, stacklevel=2)
            
        # FINALLY, BOUND dtime in case of edge cases / escapes (and set warnings)
        if dt_RESULT <= dt_min:
            dt_RESULT = dt_min
            warnings.warn("dt below minimum reasonable limit, reset to {} nanosec".format(
                round(dt_min*10**6)), category=RuntimeWarning, stacklevel=2)
        elif dt_RESULT >= dt_max:
            dt_RESULT = dt_max
            warnings.warn("dt above maximum reasonable limit, reset to {} nanosec".format(
                round(dt_max*10**6)), category=RuntimeWarning, stacklevel=2)
                
        return dt_RESULT    
       
    
    def iterate(self, ConvectionWatts, dtime:float=10**-7, 
                dynamic_dt=False, debug=False):
        """
        Iterate a Droplet:
         * INPUT Droplet Parameters (Listed Above)
         * INPUT dtime (turns Watts into Joules)
         * INPUT External Heat (Watt or Joules) from Convection?
         * OPTION Dynamic dtime OVERRIDE to input dtime (NEED TO MAKE IT CLEAR
                                                         WHEN USED!!!)
        
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
        
        # CHECK and/or REASSIGN delta Time value for Stability
        if dynamic_dt:
            # Simply calculate stable dtime and use the 0.5 margin for safety
            dtime = self._stable_dtime(margin=0.5)
            warnings.warn("NOTE: RE-ASSIGNING dtime inside of Model. DO NOT TRUST External Time Counters!"+
                          "\n\t: Use SprayDroplet.Age for Trusted Time Tracking!",
                            category=RuntimeWarning, stacklevel=2)
        else:
            self._stable_dtime(dt_check = dtime)
        
        # Calculate Heat Flux
        self.dQdt_Flux = self.k_upper[:-1] * (self.T_degC[1:] - self.T_degC[:-1])
        #       FLUX for an area is k_area (k*A/dr) * DeltaT
        #                        is for each area ABOVE the datapoint (K Upper). 
        #                        (No area below point zero, don't care about area above final point)
        #       Therefore flux will be ADDED to the one below (No change at top)
        #       And flux will be SUBTRACTED from the one at top (No change at bottom)
        #       Then we tack on the external convection wattage
        # NOTE FLUX IS LENGTH N-1 THEREFORE CANT BE NAIIVELY ADDED TO RECORDS!
        self.dQdt_net        = np.zeros_like(self.T_degC)
        self.dQdt_net[:-1]  += self.dQdt_Flux       #ADD flux to all but top
        self.dQdt_net[1:]   -= self.dQdt_Flux       #LOSE flux to all but bottom
        self.dQdt_net[-1]   += ConvectionWatts      #ADD Convective to top
        
        self.Qnet_vox       += self.dQdt_net * dtime        # dJ = W * sec
        self.Qnet_drop      += ConvectionWatts * dtime      # dJ = W * sec
        self.Age            += dtime                        # sec = sec + sec
        self.last_dt        = dtime # Just record this for debugging? 
        self.last_Conv_W    = ConvectionWatts # Record convection watts 
        
        
        # BASED ON Q Net, Update TEMPERATURE!
        # Interpolate Temperature where wet based on Cp wet
        self.T_Wet = (self.Qnet_vox <= self.H_dT_vox) * (self.initTemp + self.Qnet_vox/self.Cp_voxel_wet)
        self.T_Boil = (self.Qnet_vox > self.H_dT_vox) * (self.Qnet_vox < (self.H_dT_vox + self.Hvap_vox) ) * self.vapTemp
        self.T_Dry = (self.Qnet_vox >= (self.H_dT_vox + self.Hvap_vox) ) * (self.vapTemp + (self.Qnet_vox-(self.H_dT_vox + self.Hvap_vox))/self.Cp_voxel_dry)
        self.T_degC = self.T_Wet + self.T_Boil + self.T_Dry
        
        # IF this logic works, use same boolean math (STILL ON ENERGY) to Reassign conductivity K values?
        # Clean up later - Now using ONLY K UPPER.
        # HOWEVER. IN CASE of T_Boil, INTERPLOLATE K between Wet and Dry states??? 
        # DECISION: for now, just switch once FULLY DRY)
                #self.k_lower = ((self.T_degC > self.vapTemp) * (self.a_lower/self.dR * self.dry_condK)
                #               +(self.T_degC <= self.vapTemp) * (self.a_lower/self.dR * self.wet_condK))
        self.k_upper = ((self.Qnet_vox >= (self.H_dT_vox + self.Hvap_vox) ) * (self.a_upper/self.dR * self.dry_condK)
                       +(self.Qnet_vox < (self.H_dT_vox + self.Hvap_vox) ) * (self.a_upper/self.dR * self.wet_condK))
        
        # Calculate Report Moisture Content based on Q Net
        self.m_liquid = ((self.Qnet_vox <= self.H_dT_vox) * self.m_liquid_init +
                         (self.Qnet_vox > self.H_dT_vox) * (self.Qnet_vox <=(self.H_dT_vox + self.Hvap_vox))*(
                                self.m_liquid_init - (self.Qnet_vox - self.H_dT_vox)/self.vapQ) +
                         (self.Qnet_vox >(self.H_dT_vox + self.Hvap_vox))* 0
                        )
        self.m_voxel = self.m_solid+self.m_liquid
        self.solidsPct = self.m_solid / self.m_voxel
        
        
        self.T_degC_avg = np.average(self.T_degC, weights=self.m_voxel)
        self.solidsPct_avg = np.average(self.solidsPct, weights = self.v_voxel)
        self.isBurnt = bool(max(self.T_degC) >= self.dry_damageTemp)
        self.isHappy = bool(min(self.T_degC) >= self.dry_targetTemp)
    
    def _init_history(self, n_records):
        """
        Initialize the droplet's self.history object arrays... 
        INPUT: 
            n_timesteps:    float, number of columns to make for timesteps. 
                            (if too few, the write_history function will be slower?)
        Returns
        -------
        self.history:       DropletHistory class object, with filled out R array
                            and the correct 2D and 1D np.NaN arrays.
                            
                            
        TdegC_profile: np.ndarray   # 2D radial Temperature vs time
        Solid_profile: np.ndarray   # 2D radial Solids %Pct vs time
        TdegC_avg: np.ndarray       # 1D History of Average Temperature
        Solid_avg: np.ndarray       # 1D History of solids Pct average
        # ^Removed the "Surface" and "Midpoint" data here... 
        #       It's just the [0] and [-1] radius levels in the existing arrays. 
        #       A future plotter may also track the 0.25 and 0.75R points, don't make unneeded ones
        Flux_profile: np.ndarray    # 2D radial Heat flux (Joules) IN from UPPER Surface
        dQdt_profile: np.ndarray    # 2D radial Net Heat change per time step
        Qnet_profile: np.ndarray    # 2D radial Net Heat Accounting (Integral of dQdt)
        # droplet history summary data points of interest
        n_records: int
        dt_min_used: int                     
                            
                            
        """
        n_vox = len(self.r_vals)
        self.history.r              = self.r_vals                  # 1D plotting Axis: radius of voxels
        self.history.t              = np.full(n_records, np.nan)   # 1D plotting Axis: record of timesteps
        
        self.history.TdegC_profile  = np.full((n_records, n_vox), np.nan)   # 2D radial Temperature vs time
        self.history.Solid_profile  = np.full((n_records, n_vox), np.nan)   # 2D radial Solids %Pct vs time
        self.history.TdegC_avg      = np.full(n_records, np.nan)   # 1D History of Average Temperature
        self.history.Solid_avg      = np.full(n_records, np.nan)   # 1D History of Solids Pct
               
        self.history.Flux_profile   = np.full((n_records, n_vox), np.nan)    # 2D radial Heat flux (Joules) IN from Upper Surface
        self.history.dQdt_profile   = np.full((n_records, n_vox), np.nan)    # 2D radial Net Heat change per time step
        self.history.Qnet_profile   = np.full((n_records, n_vox), np.nan)    # 2D radial Net Heat Accounting (Integral of dQdt)
        
        self.history.dt_min_used    = np.full(n_records, np.nan)   # history of solver - recording dt timesteps used in simulations
    
    def _update_history(self, time_index=np.nan):
        """
        Record a row of results into the self.history object
        
        In a looping model, each timestep should be recorded into the NEXT time index location. 
        
        If time_index = nan, then search for the first index in time (history.t) that is np.NaN
                               If there is a time[i] = np.NaN, that i is the index. 
                               Othewise APPEND all with a new row and give a warning about that... 
                              
        For now, if we need to APPEND, I'm using np.append and or np.vstack'

        """
        if np.isnan(time_index):
            #Find index i so that we overwrite the next timestamp
            i=0
            for index in self.history.t:
                if np.isnan(index):
                    break  # Stop searching, i is the index. 
                else:
                    i += 1 # Advance i to next. (if NO np.NAN, i will end at len+1) 
                    
        # Note: FLux Record is length N-1 (because it's per between-space area)
        #   So for the RECORD, append the convection heat to the outside!
        fluxRecord = np.append(self.dQdt_Flux, self.last_Conv_W)             
                    
                    
        if i <= len(self.history.t):
        # ^ There is probably a better way to do all this, and I'm not sure that just adding a new values at a missing index 
            
            self.history.t[i]               = self.Age   # 1D plotting Axis: record of timesteps
    
            self.history.TdegC_profile[i,:] = self.T_degC        # 2D radial Temperature vs time
            self.history.Solid_profile[i,:] = self.solidsPct     # 2D radial Solids %Pct vs time
            self.history.TdegC_avg[i]       = self.T_degC_avg    # 1D History of Average Temperature
            self.history.Solid_avg[i]       = self.solidsPct_avg # 1D History of Solids Pct
            
            self.history.Flux_profile[i,:]  = self.fluxRecord     # 2D radial Heat flux (Joules) IN from Upper Surface
            self.history.dQdt_profile[i,:]  = self.dQdt_net      # 2D radial Net Heat change per time step
            self.history.Qnet_profile[i,:]  = self.Qnet_vox      # 2D radial Net Heat Accounting (Integral of dQdt)
            
            self.history.Qnet_droplet[i]    = self.Qnet_drop     # 2D radial Net Heat Accounting (Integral of dQdt)
            self.history.dt_min_used[i]     = self.last_dt       # 2D history of solver - recording dt timesteps used in simulations
            
        else: # NEED TO APPEND
            warnings.warn("APPENDING RECORDS - Model may be Slower Due to Lack of Allocated Rows?",
              category=RuntimeWarning, stacklevel=2)
            self.history.t   = np.append(self.history.t, self.Age, axis=0)   # 1D plotting Axis: record of timesteps
    
            self.history.TdegC_profile = np.vstack(self.history.TdegC_profile, self.T_degC, axis=0)        # 2D radial Temperature vs time
            self.history.Solid_profile = np.vstack(self.history.Solid_profile, self.solidsPct, axis=0)     # 2D radial Solids %Pct vs time
            self.history.TdegC_avg     = np.append(self.history.TdegC_avg, self.T_degC_avg)            # 1D History of Average Temperature
            self.history.Solid_avg     = np.append(self.history.Solid_avg, self.solidsPct_avg)         # 1D History of Solids Pct
            
            self.history.Flux_profile  = np.vstack(self.history.Flux_profile, fluxRecord, axis=0)         # 2D radial Heat flux (Joules) IN from Upper Surface
            self.history.dQdt_profile  = np.vstack(self.history.dQdt_profile, self.dQdt_net, axis=0)      # 2D radial Net Heat change per time step
            self.history.Qnet_profile  = np.vstack(self.history.Qnet_profile, self.Qnet_vox, axis=0)      # 2D radial Net Heat Accounting (Integral of dQdt)
            
            self.history.Qnet_droplet  = np.append(self.history.Qnet_droplet, self.Qnet_drop, axis=0)  # 1D history of droplet total Q
            self.history.dt_min_used   = np.append(self.history.dt_min_used, self.last_dt, axis=0)     # 1D history of solver - recording dt timesteps used in simulations
        
        
    def run_simulation(self, targetAge:float=1.0, deltaTime:float=50*10**-9, 
                        dynamicTime=False, logInterval:float=0.025,
                        otherLogging=False):
        """
        Adding a simplified Model Run and Data Collection scheme in order to
        capture 
        INPUT:
            targetAge:  float, iterate the model until self.Age >= targetAge
            deltaTime:  float, dt step for each chunk of time solving
            dynamicTime:BOOL, whether or not to re-calculate time step based
                            on the Euler and Fourier math we have done so far.
            logInterval:float, seconds time-step on which to log a new output row
            otherLogging:BOOL, Hard coded recommendations for logging other
                            points of interest such as First/Last Boil, Dry, 
                            TargetTemp, and DamageTemp, etc. 
        
        NOTE: NEEDS A SPRAY DRYER FOR THE INPUT AND/OR AIRFLOW
                TO DETERMINE THE HEAT PROFILE
        
        ACTIONS: 
            
        OUTPUT:
            runResult:  DropletHistory Class (see above?)
                        1D labels for Radius and Time
         
         
        """
        print("Empty function for now - need Sprayer Input to proceed")
        

        
# Generate Defaults

        
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
 
 #---------------------------------------------------------------------------
 #-NOTE ON MINIMUM D-TIME THEORY---------------------------------------------
 #-----"BIOT NUMBER" - Ratio of Convective vs Radial Heat Cond.: Bi = L*h/k--
 #-------------------- (For Shell, L ~ dR, so inherently smaller by dR/R)----
 #---------------------------------------------------------------------------
 #-----"FOURIER ##" - Would "Real World" heat be outpacing the time step? ---
 #------------------- AKA "NO JUMPING OVER VOXELS!"... ----------------------
 #--------------------- Fo = alpha * dt / dR^2-(for spherical coordinates)---
 #----------------------@{alpha} = k / rho*Cp = [Wm/m2K] / [g/m3][J/gK]------
 #------------------------------ = [m2/sec]----------------------------------
 #--------------------dt must stay below ~dR2/Alpha (Thus keeping Fo < 1.0!)-
 #------------------------------ OR geometry rules call for <0.5 or <0.333...  
 #---------------------------------------------------------------------------
 #---------------------------------------------------------------------------
 #---------------------------------------------------------------------------
 

"""










