#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:06:25 2020

@author: sotzee
"""
import numpy as np
import scipy.optimize as opt
from joblib import Parallel, delayed
import unitconvert
import eos_quarkyonic
import eos_potential
from config import saturation_density,u_max,N1,N2,path,args,eos_shape,L,BE,Sv,gamma,num_cores,nB_grid,args_Lattimer_pnm_init,s_k
import os
def ensure_dir(path,dir_name):
    try:
        os.stat(path+dir_name)
    except:
        os.mkdir(path+dir_name)

# Define PNM potential fixed by: BE+Sv, L, gamma
Potential_Lattimer_PNM_list=[]
for L_i in L:
    args_Lattimer=[]
    args_Lattimer_pnm=np.concatenate((opt.root(eos_potential.fit_lattimer_pnm,args_Lattimer_pnm_init,args=([unitconvert.m_n_MeV+BE+Sv,L_i,gamma],),method='hybr').x,[gamma]))
    Potential_Lattimer_pnm=eos_potential.Potential_single(args_Lattimer_pnm,eos_potential.syms_Lattimer,eos_potential.V_Lattimer_expr)
    Potential_Lattimer_PNM_list.append(Potential_Lattimer_pnm)

# use joblib to take advantage of mutiple cpu in computer.
def main_parallel_unsave(Calculation_i,parameter_list,other_args=[],verbose=1):
    Output=Parallel(n_jobs=num_cores,verbose=verbose)(delayed(Calculation_i)(parameter_i,other_args) for parameter_i in parameter_list)
    return np.array(Output)
def Calculation_creat_EOS(eos_args_args_array,other):
    return eos_quarkyonic.EOS_Quarkyonic_Potential(eos_args_args_array,potential=other,s_k=s_k,defaut_u_max=u_max)

# Define EOSs of all parameter sets.
eos_flat=[]
for Potential_Lattimer_PNM_i,L_i in zip(Potential_Lattimer_PNM_list,L):
    print('Calculating %d*%d configuration with L=%.2f MeV'%(eos_shape[1],eos_shape[2],L_i))
    eos_flat.append(main_parallel_unsave(Calculation_creat_EOS,args[1:,0].reshape((2,-1)).transpose(),other_args=Potential_Lattimer_PNM_i))
print('EOS calculated successfully!\n')
eos_flat=np.array(eos_flat).flatten()
eos=eos_flat.reshape(eos_shape)
# eos contains 3D array of EOS of shape (N_L, N_Lambda, N_nt).
# eos_flat is 1D array of all EOS of lenth N_L*N_Lambda*N_nt.
# e.g. if args = np.mgrid[20:80:4j,300:2500:23j,0.2:0.6:21j],
# eos[1,5,5]=eos_flat[1*23*21+5*21+5], for L=40 MeV, Lambda=800 MeV, nt=0.3 fm^-3


#generate a grand EOS table with regular k_Fn grid.(Analytical)
print('\nGenerating a grand EOS table with regular k_Fn grid...')
table_dir_name='EOS_table_regular_kFn'
ensure_dir(path,table_dir_name)
eos_crust=eos_flat[0].crust_eos_array
eos_crust=np.concatenate((eos_crust,eos_flat[0].eosCs2(eos_crust[2])[np.newaxis]),axis=0)
eos_crust_header='baryon number density nB (fm-3), \t energy density epsilon (MeV fm-3), \t pressure p (MeV fm-3), \t sound speed square cs2/c2'
np.savetxt(table_dir_name+'/crust_eos.txt',eos_crust.transpose(),fmt='%.8g',header=eos_crust_header)
kFn_n_e_p_cs2=np.zeros((5,N1+N2-1,np.prod(eos_shape)))
k0_kt_kmax=np.zeros((3,np.prod(eos_shape)))
for i in range(np.prod(eos_shape)):
    kFn_n_e_p_cs2[0,:,i]=eos_flat[i].k_Fn_array
    kFn_n_e_p_cs2[1:4,:,i]=eos_flat[i].eos_array_high
    kFn_n_e_p_cs2[4,:,i]=eos_flat[i].eosCs2(kFn_n_e_p_cs2[3,:,i])
    k0_kt_kmax[:,i]=np.array([eos_flat[i].k_F0,eos_flat[i].k_Ft,eos_flat[i].k_Fmax])
kFn_n_e_p_cs2_name=['neutron_fermi_momentum','baryon_number_density','energy_density','pressure','sound_speed_square']
for i in range(len(kFn_n_e_p_cs2_name)):
    np.savetxt(table_dir_name+'/core_'+kFn_n_e_p_cs2_name[i]+'.txt',kFn_n_e_p_cs2[i].transpose(),fmt='%.8g')
k0_kt_kmax_header='neutron upper Fermi momentums (MeV) at: core-crust transition, quark drip transition, and 12*saturation density.'
np.savetxt(table_dir_name+'/k0_kt_kmax.txt',k0_kt_kmax.transpose(),fmt='%.8g',header=k0_kt_kmax_header)
print('Table saved in dir \'./'+table_dir_name+'\'')


#generate a grand EOS table with regular baryon number density grid.(Intepolation)
print('\nGenerating a grand EOS table with regular baryon number density grid...')
table_dir_name='EOS_table_regular_nB'
ensure_dir(path,table_dir_name)
eos_crust=eos_flat[0].crust_eos_array
eos_crust=np.concatenate((eos_crust,eos_flat[0].eosCs2(eos_crust[2])[np.newaxis]),axis=0)
eos_crust_header='baryon number density nB (fm-3), \t energy density epsilon (MeV fm-3), \t pressure p (MeV fm-3), \t sound speed square cs2/c2'
np.savetxt(table_dir_name+'/crust_eos.txt',eos_crust.transpose(),fmt='%.8g',header=eos_crust_header)
e_p_cs2=np.zeros((3,len(nB_grid),np.prod(eos_shape)))
for i in range(np.prod(eos_shape)):
    e_p_cs2[1,:,i]=eos_flat[i].eosPressure_frombaryon(nB_grid)
    e_p_cs2[0,:,i]=eos_flat[i].eosDensity(e_p_cs2[1,:,i])
    e_p_cs2[2,:,i]=eos_flat[i].eosCs2(e_p_cs2[1,:,i])
e_p_cs2_name=['energy_density','pressure','sound_speed_square']
for i in range(len(e_p_cs2)):
    np.savetxt(table_dir_name+'/core_'+e_p_cs2_name[i]+'.txt',e_p_cs2[i].transpose(),fmt='%.8g')
np.savetxt(table_dir_name+'/core_nB_grid.txt',nB_grid,fmt='%.8g')
print('Table saved in dir \'./'+table_dir_name+'\'')



# e.g. if args = np.mgrid[20:80:4j,300:2500:23j,0.2:0.6:21j],
# eos[1,5,5]=eos_flat[1*23*21+5*21+5], for L=40 MeV, Lambda=800 MeV, nt=0.3 fm^-3
eos_to_check=eos[1,5,5]
# To get [n, epsilon, p] array of a EOS, modify N1, N2 to vary the size of the array
eos_array_test = eos_to_check.eos_array
n_array=eos_array_test[0]
epsilon_array=eos_array_test[1]
p_array=eos_array_test[2]

# Analyical solution for PNM quarkyonic matter is aviable given k_Fn as prime
# argument. If one need e.g. pressure as a function of baryon density,
# intepolation of array is needed. EOS has been intepolated already.
print('\nHere is the example: eos[1,5,5]')
ns=saturation_density
ps=eos_to_check.eosPressure_frombaryon(ns)    #get pressure in MeV/fm^-3
es=eos_to_check.eosDensity(ps)                #get energy density in MeV/fm^-3
cs2=eos_to_check.eosCs2(ps)                   #get sound speed square
chempo=eos_to_check.eosChempo(ps)             #get chemical potential
#ns=eos_to_check.eosBaryonDensity(ps)
print('at nB=%.6f fm-3f, e=%.3f MeV fm-3, p=%.6f, cs^2=%.6f'%(ns,es,ps,cs2))
print('L=%.6f MeV'%(3*ps/ns))
print('BE+Sv=%.6f MeV'%(es/ns-unitconvert.m_n_MeV))