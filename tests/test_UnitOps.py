# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 15:51:48 2026

@author: dhooke
"""
import numpy as np
import pytest

from spraypy import InletAir, SlurryData, SprayDryer, SprayDroplet


#------------------------------------------------------------------------------
#---SECTION: DROPLET GEOMETRY--------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def test_outer_edge_Radius_OffByOne():
    dryer = SprayDryer(dropR = 100*10**-6) #simple Spray dryer with 100um dropR
    droplet = SprayDroplet(SlurryData(), dryer, dR = 5*10*-6) #5um voxels
    outer_edge = droplet.r_vals[-1] # + droplet.dR <-- Not Off-by-one error? 
    assert outer_edge == pytest.approx(dryer.dropR)
    
def test_droplet_volume_mass_correct():
    dryer = SprayDryer(dropR = 100*10**-6) #simple Spray dryer with 100um dropR
    droplet = SprayDroplet(SlurryData(),dryer, dR = 5*10*-6) #5um voxels
    volume_sum = np.sum(droplet.v_voxel)
    #mass_sum = np.sum(droplet.m_voxel)
    #assert volume_sum == pytest.approx(4/3*np.pi*dryer.dropR**3)
    assert volume_sum == pytest.approx(4/3*np.pi*(dryer.dropR+droplet.dR)**3)
    
    #assert mass_sum   == pytest.approx(droplet.m_liquid + droplet.m_solid)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------