# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 18:03:02 2016

@author: sotzee
"""

from astropy.constants import M_sun
import scipy.constants as const

c=100*const.c
G=1000*const.G
e=1e7*const.e
hbar=const.hbar*1e7
mass_sun=1000*M_sun.value
mass_per_baryon=const.m_n*1000
m_p_MeV=const.m_p*const.c**2/(1e6*const.e)
m_n_MeV=const.m_n*const.c**2/(1e6*const.e)
unitMeVfm=((1e6*e/hbar/c)/1e13)**3
unitPressure=(1e6*e)**4/(hbar*c)**3
unitDensity=(1e6*e)**4/(hbar*c)**3/c**2
unitBaryonDensity=1e-39
unitMass=mass_sun
Gc2=G/c**2

def toPressure(pressure_before,unit_before):
    if(unit_before=='mevfm' or unit_before=='mevfm3' or unit_before=='mevfm-3'):
        return pressure_before/unitMeVfm*unitPressure
    if(unit_before=='mev4'):
        return pressure_before*unitPressure
    if(unit_before=='mev'):
        return pressure_before**4.0*unitPressure
    if(unit_before=='fm-4'):
        return pressure_before/unitMeVfm**(4./3)*unitPressure
        
def toDensity(density_before,unit_before):
    if(unit_before=='mevfm' or unit_before=='mevfm3' or unit_before=='mevfm-3'):
        return density_before/unitMeVfm*unitDensity
    if(unit_before=='mev4'):
        return density_before*unitDensity
    if(unit_before=='mev'):
        return density_before**4.0*unitDensity
    if(unit_before=='fm-4'):
        return density_before/unitMeVfm**(4./3)*unitDensity

def toBaryonDensity(density_before,unit_before):
    if(unit_before=='mevfm' or unit_before=='mevfm3' or unit_before=='mevfm-3'):
        return density_before/unitBaryonDensity
        
def toMevfm(before,unit_before):
    if(unit_before=='pressure'):
        return before*unitMeVfm/unitPressure
    if(unit_before=='density'):
        return before*unitMeVfm/unitDensity
    if(unit_before=='baryondensity'):
        return before*unitBaryonDensity
    if(unit_before=='mev4'):
        return before*unitMeVfm
        

def toMev4(before,unit_before):
    if(unit_before=='pressure'):
        return before/unitPressure
    if(unit_before=='density'):
        return before/unitDensity
    if(unit_before=='mevfm' or unit_before=='mevfm3' or unit_before=='mevfm-3'):
        return before/unitMeVfm
        
def toMev(before,unit_before):
    if(unit_before=='pressure'):
        return (before/unitPressure)**0.25
    if(unit_before=='density'):
        return (before/unitDensity)**0.25

def toRadius(before,unit_before):
    if(unit_before=='beta'):
        return mass_sun*Gc2*before/100000

def toBeta(before,unit_before):
    if(unit_before=='mass'):
        return mass_sun*Gc2*before/100000
    if(unit_before=='radius'):
        return mass_sun*Gc2*before/100000
        
#test:
#print('test:')
#print(toMevfm(toPressure(10,'mevfm'),'pressure'))
#print(toMevfm(toDensity(10,'mevfm'),'density'))
#print(toMev4(toPressure(10,'mev4'),'pressure'))
#print(toMev4(toDensity(10,'mev4'),'density'))
#print(toMev(toPressure(10,'mev'),'pressure'))
#print(toMev(toDensity(10,'mev'),'density'))