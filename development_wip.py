# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:14:28 2026

@author: dhooke
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close("all") #Get rid of plots from last run of program
#import scipy as sp
from dataclasses import dat

#Use this to report better values
from math import log10, floor
def round_2(x): return round(x, 2-int(floor(log10(abs(x)))))




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
dt=0.001       # dt value, in seconds (0.1 milisecond)
dr=10*10**(-6) # dr value, in meters (10 microns, 100 points)
dtime=dr*10      # dt value, in seconds (0.1 milisecond) Not used?
t_vals=[] #New Idea. Fix delta T (at the surface), and simply RECORD the time of each step!

Plots_To_Output_T=[] #Initialize list of models
Titles_To_Output_T=[]
Rs_To_Output_T=[]
tcount_Out=[]
for Ti in T_Settings:
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
