# -*- coding: utf-8 -*-
"""
Spray.py Codebase

Spray Dryer Simulation and Modeling in Python
2020 U.Washington Masters ChemEMath Project, David Hurt Hooke

Initial Jupyter notebook focuses on analytical mathematics of 
    spherical heat transfer in an arbitrary small sphere of given 
    thermal properties (to initially simulate water and slurry)
Boundary conditions of (at first) constant external temperature 
    and heat transfer coefficient represent arbitrarilly large 
    Spray Dryer and/or minimum flow rate of particles. 
Initial results show and compare the influence of size and heat 
    transfer properties on the uniformity of heating (for sensitive
    materials such as food or medicine, it is important to dry 
    the middle without over-heating the outside material!)
    
The second stage added the heat of vaporization, where the droplet
    actually reaches a boiling point (constant for now?) and therefore
    packets will continue to transfer heat down any thermal gradient but 
    will not keep heating up until they have boiled, and then they 
    have a new thermal mass and conductivity of the solid(and gas) phase.
    (This does not deal with any gas-transport kinetics? hmm).

The third stage , the final one considered but not finished in that 
    project was to return  to the spray dryer itself, where a particular 
    feed rate and setup will impact the atmosphere an average droplet 
    is exposed to during its existance, both by the total heat injected
    into the system by oncoming air and by the momentum influence on 
    Residence Time in the dryer... (This will need fudgig as flow is 
                                    VERY complicated and swirley! )

Other things to consider for future work: 
    * Building a library of materials / recipes to spray
    * Learning about pneumatic dust conveyance for particle handling
    * Other changes inside particles? (Hollow sphere results?)

@author: dhooke
Converting from Jupyter Notebook in Fall of 2026 as a fun project
to recall old times, and explore publishing engineering-toolkits online
using StreamLit App service.

If this ever ends up being useful, by all means someone with more time 
and expertise should take this inspiration and actually do it right!

"""
#----------------------------------------------------


"""
------ ORGANIZATION + PLANNING -----------------
For the moment I'm copying directly from my old jupyter work
to get something that runs, and we can work on turning it into 
better Pythonic functions later. 

I still don't really understand "Classes" in python, but is this a thing
where we would make a "Spray Dryer" class to set up the macro physics, 
a "Recipe" class to set up the chemistry of the droplet, and then a "Droplet" 
class to take the recipe and put it into the spray dryer? I totally don't know.'

In any case, I'll function-alize what I can and convert the existing into 
tests to run it.... 
"""




"""
-------------------------------------------------
------- INITIAL JUPYTER BELOW HERE --------------
-------------------------------------------------
------- Comment out when Obsolete ---------------
-------------------------------------------------
"""
"""
------A-------
-----A-A------
----A---A-----
---A--A--A----
--A-------A---
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close("all") #Get rid of plots from last run of program
import scipy as sp

#Use this to report better values
from math import log10, floor
def round_2(x): return round(x, 2-int(floor(log10(abs(x)))))

#Describe the Droplet
R=300*10**(-6) # Droplet Radius, meters (eg. 500 microns)
h=1500          # Watts/m^2 K (convective heat transfer. a "Typical" value for flowing hot air)
k=0.500         # Watts/m K (Thermal Conductivity (Assumed Constant Here, looked up normal value for water...))
Cp=4.18         # J/g K     (Again, lookup value for water)
rho=1*10**(6)   # g/m^3      (Water by definition)
HVap=2256.4     # J/g to Evaporate water. NOTE: is 500x the Q for 1K dT - more than 5x the Q for 100C change

#Describe the System
Ti=300 + 273.16 #K, heat of system (Maybe make f(t) eventually)
T0=25  + 273.16 #K, Initial Temperature of water (Constant at all points)

time=0.5 # Seconds in the spray dryer (max time to solve for)
dt=0.001       # dt value, in seconds (0.1 milisecond)
t_vals=[] #New Idea. Fix delta T (at the surface), and simply RECORD the time of each step!

dr=10*10**(-6) # dr value, in meters (10 microns, 100 points)
r_vals=np.arange(0,R+dr,dr)
A_Vals=4*np.pi*r_vals**2 #Area of the shell at r 
V_Vals=r_vals*0
for i,r in enumerate(r_vals[1:]):
    V_Vals[i+1] = 4/3*np.pi*r**3 - 4/3*np.pi*r_vals[i]**3 #Volume of the shell at r - 
H_Vals=rho*V_Vals*HVap


Temps=r_vals*0+T0 #Set initial Temperature profile = T0 is constant?
dTdt=r_vals*0 #Set dTdt 
Biot=(R/3)*h/k

# fig, ax = plt.subplots(figsize=[10,5])
# ax.set_ylim([0,300])
# ax.set_ylabel("Temperature *C")
# ax.set_xlabel("Radius (Microns)")
# ax.set_title("Heating with T="+str(round(Ti-273.16))+'*C   Biot='+str(round_2(Biot)))
# ax.axhspan(125, 175, alpha=0.1, color='black')

Plot_T=np.linspace(273,Ti+10,50) #Plot n solutions
azm=np.linspace(0,2*np.pi,365)
r_V,Th_V=np.meshgrid(r_vals,azm)
Plots,a=np.meshgrid(r_vals,Plot_T)
#for t in t_vals:
#     dt=0.0001
#     #For each timestamp, first do the balance on r at the edge. 
    
dT=0.001 #Kelvin per step for initial solution
p=0 # Index of temperatures to plot.
t=10**-15 # Initiate Time-Counter

T_S= Temps[-1]
while T_S <Ti*0.90:
    
    T_S= Temps[-1] #Record surface temperature, will decide which to plot
    
    if T_S > Plot_T[p]:
        Plots[p,:]=Temps
        p=p+1
        #print(str(round_2(t)*10**3)+" milliSeconds")
    
    for i,r in enumerate(r_vals): #Loop through all r values to find dT/dt
        if r==0:
            Q_out=0         #Boundary condition at middle - no place for that temperature to go.
        else:
            Q_out=(4*np.pi*(r)**2)*k*(Temps[i]-Temps[i-1])/dr    #Conduction through shell at r

        if r==max(r_vals):
            Q_in = (4*np.pi*(r**2))*h*(Ti-Temps[i])    #Convection at edge of shell
        else:
            Q_in = (4*np.pi*((r+dr)**2))*k*(Temps[i+1]-Temps[i])/dr #Conduciton through shell at -r

        Q = Q_in - Q_out

        V=4/3*np.pi*((r+dr)**3-(r)**3)

        dTdt[i]= (Q) / (rho*V*Cp) # K/s,   (Watts) / (g/m3 * m3 * J/gK)
        
    if dTdt[-1]>0:
        deltaT=dT* dTdt/dTdt[-1] #CALCULATE the time step as the time required for a dT-sized change at the surface
        t=t+dT/dTdt[-1] #Capture the time elapsed
    #small delta T was required for smooth startup. I think we can use faster dT steps once we're in the weeds!
    if T_S>300:
        dT=0.05
    elif T_S>400:
        dT=0.1
    
    Temps=Temps+deltaT

fig, ax = plt.subplots(figsize=[10,5])
ax.set_ylim([0,300])
ax.set_ylabel("Temperature *C")
ax.set_xlabel("Radius (Microns)")
ax.set_title("Heating with T="+str(round(Ti-273.16))+'*C   Biot='+str(round_2(Biot)))
ax.axhspan(125, 175, alpha=0.1, color='black')
for i,T in enumerate(Plot_T):    
    ax.plot(r_vals*10**6,Plots[i,:]-273.15)



"""
----BBBB---------
---BBB--BB-------
---BBB---BB------
---BB-BBBB-------
---BBB---BB------
---BBBBBB--------
"""

# DIFFERENT DROPLET SIZES 


#Describe the Droplet
R=300*10**(-6) # Droplet Radius, meters (eg. 500 microns)
R_Sizes=[100*10**(-6),300*10**(-6),500*10**(-6),900*10**(-6)]

h=1000          # Watts/m^2 K (convective heat transfer. a "Typical" value for flowing hot air)
k=0.500         # Watts/m K (Thermal Conductivity (Assumed Constant Here, looked up normal value for water...))
Cp=4.18         # J/g K     (Again, lookup value for water)
rho=1*10**(6)   # g/m^3      (Water by definition)
HVap=2256.4     # J/g to Evaporate water. NOTE: is 500x the Q for 1K dT - more than 5x the Q for 100C change

#Describe the System
Ti=300 + 273.16 #K, heat of system (Maybe make f(t) eventually)
T0=25  + 273.16 #K, Initial Temperature of water (Constant at all points)

time=0.5 # Seconds in the spray dryer (max time to solve for)
dt=0.0001       # dt value, in seconds (0.1 milisecond)
t_vals=[] #New Idea. Fix delta T (at the surface), and simply RECORD the time of each step!

Plots_to_Make_R=[] #Initialize list of "Plots" to make
Rs_to_Plot=[]
Titles=[]
for R in R_Sizes:
    dr=10*10**(-6) # dr value, in meters (10 microns, 100 points)
    r_vals=np.arange(0,R+dr,dr)
    A_Vals=4*np.pi*r_vals**2 #Area of the shell at r 
    V_Vals=r_vals*0
    for i,r in enumerate(r_vals[1:]):
        V_Vals[i+1] = 4/3*np.pi*r**3 - 4/3*np.pi*r_vals[i]**3 #Volume of the shell at r - 
    H_Vals=rho*V_Vals*HVap


    Temps=r_vals*0+T0 #Set initial Temperature profile = T0 is constant?
    dTdt=r_vals*0 #Set dTdt 
    Biot=(R/3)*h/k


    Plot_T=np.linspace(273,Ti+10,30) #Plot n solutions
    #Prepare things for Plotting
    azm=np.linspace(0,2*np.pi,365)
    r_V,Th_V=np.meshgrid(r_vals,azm)
    Plots,a=np.meshgrid(r_vals,Plot_T)

    #for t in t_vals:
    #     dt=0.0001
    #     #For each timestamp, first do the balance on r at the edge. 

    dT=0.0005 #Kelvin per step for initial solution
    p=0 # Index of temperatures to plot.
    t=10**-15 # Initiate Time-Counter

    T_S= Temps[-1]
    while T_S <Ti*0.95 and Temps[1]<125+273.16:
        #while Temps[1]<100+273.15: #Do until innermost has reached temperature
        T_S= Temps[-1]

        if T_S > Plot_T[p]:
            Plots[p,:]=Temps
            p=p+1
            
            #print(str(round_2(t)*10**3)+" milliSeconds")

        for i,r in enumerate(r_vals): #Loop through all r values to find dT/dt
            if r==0:
                Q_out=0         #Boundary condition at middle - no place for that temperature to go.
            else:
                Q_out=(4*np.pi*(r)**2)*k*(Temps[i]-Temps[i-1])/dr    #Conduction through shell at r

            if r==max(r_vals):
                Q_in = (4*np.pi*(r**2))*h*(Ti-Temps[i])    #Convection at edge of shell
            else:
                Q_in = (4*np.pi*((r+dr)**2))*k*(Temps[i+1]-Temps[i])/dr #Conduciton through shell at -r

            Q = Q_in - Q_out

            V=4/3*np.pi*((r+dr)**3-(r)**3)

            dTdt[i]= (Q) / (rho*V*Cp) # K/s,   (Watts) / (g/m3 * m3 * J/gK)

        if dTdt[-1]>0:
            deltaT=dT* dTdt/dTdt[-1] #CALCULATE the time step as the time required for a dT-sized change at the surface
            t=t+dT/dTdt[-1] #Capture the time elapsed
        #small delta T was required for smooth startup. I think we can use faster dT steps once we're in the weeds!
        if T_S>300:
            dT=0.01
        elif T_S>400:
            dT=0.02

        Temps=Temps+deltaT
    Plots[p,:]=Temps #Capture final point
    Titles.append("Heating with T="+str(round(Ti-273.16))+'*C   Biot='+str(round_2(Biot)) + '   TIME = ' + str(round_2(t))+" Seconds")
    Plots_to_Make_R.append(Plots) #Save this plot for later!
    Rs_to_Plot.append(r_vals)
    print("Done R={}um \t Tmin=".format(round(R*10**6))+str(round(Temps[0]-273))+ "*C    Tmax="+str(round(Temps[-1]-273))+"*C")
print("Thanks!")



for j,Plots in enumerate(Plots_to_Make_R):
    r_vals=Rs_to_Plot[j]
    fig, ax = plt.subplots(figsize=[10,3])
    ax.set_ylim([0,300])
    ax.set_xlim([0,1000])
    ax.set_ylabel("Temperature *C")
    ax.set_xlabel("Radius (Microns)")
    for i,y in enumerate(Plot_T):
        ax.plot(r_vals*10**6,Plots[i,:]-273.15,label="time=" +str(round_2(t)*10**9)+"ns")
    ax.set_title(Titles[j])
    if (np.max(Plots)-273.15)>175:        
        ax.axhspan(115, 175, alpha=0.1, color='red')
    else:
        ax.axhspan(115, 175, alpha=0.1, color='green')   
"""
-----CCC------
---CC---CC----
--CC----------
--CC----------
--CCC---CC----
---CCCCC------
"""
# DIFFERENT INPUT TEMPERATURES


#Describe the Droplet
R=300*10**(-6) # Droplet Radius, meters (eg. 500 microns)
T_Settings=[150+273.16, 250+273.16, 350+273.16, 450+273.16]


h=1000          # Watts/m^2 K (convective heat transfer. a "Typical" value for flowing hot air)
k=0.500         # Watts/m K (Thermal Conductivity (Assumed Constant Here, looked up normal value for water...))
Cp=4.18         # J/g K     (Again, lookup value for water)
rho=1*10**(6)   # g/m^3      (Water by definition)
HVap=2256.4     # J/g to Evaporate water. NOTE: is 500x the Q for 1K dT - more than 5x the Q for 100C change

#Describe the System
#Ti=300 + 273.16 #K, heat of system (Maybe make f(t) eventually) 
T0=25  + 273.16 #K, Initial Temperature of water (Constant at all points)

time=0.5 # Seconds in the spray dryer (max time to solve for)
dt=0.001       # dt value, in seconds (0.1 milisecond)
t_vals=[] #New Idea. Fix delta T (at the surface), and simply RECORD the time of each step!

Plots_To_Output_T=[]
for Ti in T_Settings:
    dr=10*10**(-6) # dr value, in meters (10 microns, 100 points)
    r_vals=np.arange(0,R+dr,dr)
    A_Vals=4*np.pi*r_vals**2 #Area of the shell at r 
    V_Vals=r_vals*0
    for i,r in enumerate(r_vals[1:]):
        V_Vals[i+1] = 4/3*np.pi*r**3 - 4/3*np.pi*r_vals[i]**3 #Volume of the shell at r - 
    H_Vals=rho*V_Vals*HVap


    Temps=r_vals*0+T0 #Set initial Temperature profile = T0 is constant?
    dTdt=r_vals*0 #Set dTdt 
    Biot=(R/3)*h/k

    fig, ax = plt.subplots(figsize=[10,3])
    ax.set_ylim([0,300])
    ax.set_xlim([0,300])
    ax.set_ylabel("Temperature *C")
    ax.set_xlabel("Radius (Microns)")

    Plot_T=np.linspace(273,Ti+10,30) #Plot n solutions
    Plots,a=np.meshgrid(r_vals,Plot_T)


    #for t in t_vals:
    #     dt=0.0001
    #     #For each timestamp, first do the balance on r at the edge. 

    dT=0.005 #Kelvin per step for initial solution
    p=0 # Index of temperatures to plot.
    t=10**-15 # Initiate Time-Counter

    T_S= Temps[-1]
    while T_S <Ti*0.95 and Temps[1]<125+273.16:
        #while Temps[1]<100+273.15: #Do until 
        T_S= Temps[-1] # Surface Temperature

        if T_S > Plot_T[p]:
            p=p+1
            ax.plot(r_vals*10**6,Temps-273.15,label="time=" +str(round_2(t)*10**9)+"ns")

            #print(str(round_2(t)*10**3)+" milliSeconds")

        for i,r in enumerate(r_vals): #Loop through all r values to find dT/dt
            if r==0:
                Q_out=0         #Boundary condition at middle - no place for that temperature to go.
            else:
                Q_out=(4*np.pi*(r)**2)*k*(Temps[i]-Temps[i-1])/dr    #Conduction through shell at r

            if r==max(r_vals):
                Q_in = (4*np.pi*(r**2))*h*(Ti-Temps[i])    #Convection at edge of shell
            else:
                Q_in = (4*np.pi*((r+dr)**2))*k*(Temps[i+1]-Temps[i])/dr #Conduciton through shell at -r

            Q = Q_in - Q_out

            V=4/3*np.pi*((r+dr)**3-(r)**3)

            dTdt[i]= (Q) / (rho*V*Cp) # K/s,   (Watts) / (g/m3 * m3 * J/gK)

        if dTdt[-1]>0:
            deltaT=dT* dTdt/dTdt[-1] #CALCULATE the time step as the time required for a dT-sized change at the surface
            t=t+dT/dTdt[-1] #Capture the time elapsed
            
        #small delta T was required for smooth startup. I think we can use faster dT steps once we're in the weeds!
        if T_S>300:   # Once we pass 27 *C (Almost immediately, but still after ~100 loops)
            dT=0.01
        elif T_S>400: # Once we pass 127*C 
            dT=0.1

        Temps=Temps+deltaT
    ax.plot(r_vals*10**6,Temps-273.15,label="time=" +str(round_2(t)*10**9)+"ns")
    ax.set_title("Heating with T="+str(round(Ti-273.16))+'*C   Biot='+str(round_2(Biot)) + '   TIME = ' + str(round_2(t))+" Seconds")
    if (Temps[-1]-273.15)>175:        
        ax.axhspan(125, 175, alpha=0.1, color='red')
    else:
        ax.axhspan(125, 175, alpha=0.1, color='green')
  




"""
--DDDDD-------
--DD---DDD----
--DD-----DD---
--DD-----DD---
--DD---DD-----
--DDDDD-------
"""
# 
#Final Project for this presentation:
#Re-craft the model for constant Q steps (not t-surface)
#And add the HVap part

#Describe the Droplet
R=300*10**(-6) # Droplet Radius, meters (eg. 500 microns)
T_Settings=[150+273.16, 250+273.16, 350+273.16, 450+273.16]


h=1000          # Watts/m^2 K (convective heat transfer. a "Typical" value for flowing hot air)
k=0.500         # Watts/m K (Thermal Conductivity (Assumed Constant Here, looked up normal value for water...))
Cp=4.18         # J/g K     (Again, lookup value for water)
rho=1*10**(6)   # g/m^3      (Water by definition)
HVap=2256.4     # J/g to Evaporate water. NOTE: is 500x the Q for 1K dT - more than 5x the Q for 100C change

#Describe the System
#Ti=300 + 273.16 #K, heat of system (Maybe make f(t) eventually) 
T0=25  + 273.16 #K, Initial Temperature of water (Constant at all points)

time=0.5 # Seconds in the spray dryer (max time to solve for)
dtime=dr*10      # dt value, in seconds (0.1 milisecond) Not used?
t_vals=[] #New Idea. Fix delta T (at the surface), and simply RECORD the time of each step!

Plots_To_Output_T=[] #Initialize list of models
Titles_To_Output_T=[]
Rs_To_Output_T=[]
tcount_Out=[]
for Ti in T_Settings:
    dr=10*10**(-6) # dr value, in meters (10 microns, 100 points)
    r_vals=np.arange(0,R+dr,dr)
    A_Vals=4*np.pi*r_vals**2 #Area of the shell at r 
    V_Vals=r_vals*0
    for i,r in enumerate(r_vals[1:]):
        V_Vals[i+1] = 4/3*np.pi*r**3 - 4/3*np.pi*r_vals[i]**3 #Volume of the shell at r - 
    H_Vals=rho*V_Vals*HVap #Joules for each volume disc to evaporate 


    Temps=r_vals*0+T0 #Set initial Temperature profile = T0 is constant?
    dTdt=r_vals*0 #Set dTdt 
    Biot=(R/3)*h/k

    fig, ax = plt.subplots(figsize=[10,4])
    ax.set_ylim([0,200])
    ax.set_xlim([0,300])
    ax.set_ylabel("Temperature *C")
    ax.set_xlabel("Radius (Microns)")

    Plot_Time=np.arange(0,10,0.01) #Plot n solutions in 0.5 seconds in equal timesteps!
    Plots,a=np.meshgrid(r_vals,np.arange(0,5,dtime)) #This is how you get the data out of the loop!


    #for t in t_vals:
    #     dt=0.0001
    #     #For each timestamp, first do the balance on r at the edge. 

    dT=0.005 #Kelvin per step for initial solution
    p=0 # Index of temperatures to plot.
    tcount=0
    t=10**-15 # Initiate Time-Counter

    T_S= Temps[-1]
    while T_S <Ti*0.95 and Temps[1]<125+273.16:
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
                Q_in = (4*np.pi*(r**2))*h*(Ti-Temps[i])    #Convection at edge of shell
                
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
            t=t+dtime #Capture the time elapsed
            
        #small delta T was required for smooth startup. I think we can use faster dT steps once we're in the weeds!
#         if T_S>300:   # Once we pass 27 *C (Almost immediately, but still after ~100 loops)
#             dt=0.005
#         elif T_S>400: # Once we pass 127*C 
#             dt=0.01

        Temps=Temps+dtime*dTdt
        
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

"""
--EEEEEEEE----
--EE----------
--EE----------
--EEEEEE------
--EE----------
--EEEEEEEEE---
"""
#OK, One more step!
#Use a changing T-Inf profile! Perhaps... Tinf changes with 10% of H? Try it out. 


#Describe the Droplet
R=300*10**(-6) # Droplet Radius, meters (eg. 500 microns)
T_Settings=[200+273.16, 250+273.16, 300+273.16, 350+273.16]


h=1000          # Watts/m^2 K (convective heat transfer. a "Typical" value for flowing hot air)
k=0.500         # Watts/m K (Thermal Conductivity (Assumed Constant Here, looked up normal value for water...))
Cp=4.18         # J/g K     (Again, lookup value for water)
rho=1*10**(6)   # g/m^3      (Water by definition)
HVap=2256.4     # J/g to Evaporate water. NOTE: is 500x the Q for 1K dT - more than 5x the Q for 100C change

#Describe the System
#Ti=300 + 273.16 #K, heat of system (Maybe make f(t) eventually) 
T0=25  + 273.16 #K, Initial Temperature of water (Constant at all points)

time=0.5 # Seconds in the spray dryer (max time to solve for)
dtime=dr*10      # dt value, in seconds (0.1 milisecond) Not used?
t_vals=[] #New Idea. Fix delta T (at the surface), and simply RECORD the time of each step!

Plots_To_Output_T=[] #Initialize list of models
Titles_To_Output_T=[]
Rs_To_Output_T=[]
tcount_Out=[]
for Ti in T_Settings:
    dr=10*10**(-6) # dr value, in meters (10 microns, 100 points)
    r_vals=np.arange(0,R+dr,dr)
    A_Vals=4*np.pi*r_vals**2 #Area of the shell at r 
    V_Vals=r_vals*0
    for i,r in enumerate(r_vals[1:]):
        V_Vals[i+1] = 4/3*np.pi*r**3 - 4/3*np.pi*r_vals[i]**3 #Volume of the shell at r - 
    H_Vals=rho*V_Vals*HVap #Joules for each volume disc to evaporate 


    Temps=r_vals*0+T0 #Set initial Temperature profile = T0 is constant?
    dTdt=r_vals*0 #Set dTdt 
    Biot=(R/3)*h/k

    fig, ax = plt.subplots(figsize=[10,4])
    ax.set_ylim([0,200])
    ax.set_xlim([0,300])
    ax.set_ylabel("Temperature *C")
    ax.set_xlabel("Radius (Microns)")

    Plot_Time=np.arange(0,20,0.01) #Plot n solutions in 0.5 seconds in equal timesteps!
    Plots,a=np.meshgrid(r_vals,np.arange(0,5,dtime)) #This is how you get the data out of the loop!


    #for t in t_vals:
    #     dt=0.0001
    #     #For each timestamp, first do the balance on r at the edge. 

    dT=0.005 #Kelvin per step for initial solution
    p=0 # Index of temperatures to plot.
    tcount=0
    t=10**-15 # Initiate Time-Counter

    T_S= Temps[-1]
    Ti_Variable=Ti #Initialize changing Temperature
    while T_S <Ti*0.95 and Temps[1]<125+273.16 and tcount<len(Plots):
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
            t=t+dtime #Capture the time elapsed
            
        #small delta T was required for smooth startup. I think we can use faster dT steps once we're in the weeds!
#         if T_S>300:   # Once we pass 27 *C (Almost immediately, but still after ~100 loops)
#             dt=0.005
#         elif T_S>400: # Once we pass 127*C 
#             dt=0.01

        Temps=Temps+dtime*dTdt
        
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



"""
--------------
--------------
--------------
--------------
--------------
--------------
"""





