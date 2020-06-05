#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 21:27:36 2020

@author: sotzee
"""

import numpy as np
from multiprocessing import cpu_count

# Choose crust EOS
crust_eos_name='Sly4' #or 'Togashi'

# Saturation density
saturation_density=0.16

# Calculate EOS upto baryon number density nB=u_max*saturation_density
u_max=12

# N1 points between k0 and k_t, N2 points between k_t and k_max. 
# k_0, k_t, k_max are neutron upper Fermi momentums at
# core-crust transition, quark drip transition, and 12*saturation density.
N1=100
N2=200

# main repository
path='./'

# 3 parameters are symmetry energy slope L (MeV), shell parameter Lambda (MeV)
# and quark drip density n_t (fm^-3).
# Three numbers of each parameter means begin:end:number*1j
args = np.mgrid[20:100:5j,300:2500:23j,0.2:0.6:21j]

#L=args[0], Lambda=args[1], n_t=args[2]. eos_shape is the demension of them.
eos_shape = args[0].shape
L      = args[0,:,0,0]
Lambda = args[0,0,:,0]
n_t    = args[0,0,0,:]

# Potential is based on PNM potential in https://arxiv.org/abs/2004.08293
# SNM binding energy at saturation density
BE=-16

# Symmetry energy at saturation density
Sv=31

# Fixed PNM potential parameter gamma
gamma=5/3

# cpu_count() is to check cpu number of this machine, num_cores is the number
# of cpu to be used.
num_cores = cpu_count()-1

# baryon number density grid to generate core EOS.
nB_grid=np.linspace(0.08,1.9,183)


#Following parameters don't need to be modified.
dlnx_cs2=1e-10 #finite difference parameter in calculate speed of sound

# The initial trial for root-finding coefficients in PNM potenals
args_Lattimer_pnm_init=[-42,25]

# EOS intepolation smoothing factor s (0 for no smoothing),
# intepolation order k (2 for quardratic).
s_k=[0,2]