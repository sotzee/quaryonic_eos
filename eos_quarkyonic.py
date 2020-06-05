#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 13:38:13 2020

@author: sotzee
"""
import numpy as np
from scipy import interpolate
from scipy.optimize import root
from scipy.misc import derivative
import unitconvert
from config import crust_eos_name,saturation_density,N1,N2,path,dlnx_cs2

if(crust_eos_name=='Sly4'):
    Sly4_eos_array=np.loadtxt(path+'sly.dat',skiprows=0)
    Sly4_eos_array[:,0]=unitconvert.toMevfm(Sly4_eos_array[:,0]/1.66*1e24,'baryondensity')
    Sly4_eos_array[:,1]=unitconvert.toMevfm(Sly4_eos_array[:,1],'density')
    Sly4_eos_array[:,2]=unitconvert.toMevfm(Sly4_eos_array[:,2],'density')
    Sly4_eos_array[:,1:]=Sly4_eos_array[:,[2,1]]
    Sly4_eos_array=Sly4_eos_array.transpose()
    Sly4_match_point=[0.0622,58.99880738,0.24014206]    # Crustal EOS array is used upto this point
    crust_eos_array=Sly4_eos_array[:,Sly4_eos_array[0]<Sly4_match_point[0]*1.0001]
elif(crust_eos_name=='Togashi'):
    Togashi_eos_array=np.loadtxt(path+'togashi.dat',skiprows=4)
    Togashi_eos_array[:,0]=Togashi_eos_array[:,0]*1e-39
    Togashi_eos_array[:,1]=unitconvert.toMevfm(Togashi_eos_array[:,1],'pressure')
    Togashi_eos_array[:,2]=unitconvert.toMevfm(Togashi_eos_array[:,2],'pressure')
    Togashi_eos_array=Togashi_eos_array.transpose()
    Togashi_match_point=[0.0597043,57.00463907,0.2093527]# Crustal EOS array is used upto this point
    crust_eos_array=Togashi_eos_array[:,Togashi_eos_array[0]<Togashi_match_point[0]*1.0001]
    
class EOS_interpolation(object):
    def __init__(self,baryon_density_s,eos_array,s_k=[0,2]): #defalt s=0,k=2 equal quadratic 1d intepolation
        self.s,self.k=s_k
        n_array,energy_array,pressure_array=self.eos_array
        self.eosPressure_frombaryon = interpolate.UnivariateSpline(n_array,pressure_array, k=self.k,s=self.s)
        self.eosDensity  = interpolate.UnivariateSpline(pressure_array,energy_array, k=self.k,s=self.s)
        self.eosBaryonDensity = interpolate.UnivariateSpline(pressure_array,n_array, k=self.k,s=self.s)
    def eosChempo(self,pressure):
        return (pressure+self.eosDensity(pressure))/self.eosBaryonDensity(pressure)
    def eosCs2(self,pressure):
        return 1.0/derivative(self.eosDensity,pressure,dx=self.eosDensity(pressure)*dlnx_cs2)

class EOS_Quarkyonic_Potential(EOS_interpolation):
    baryon_density_s=saturation_density
    m_p=unitconvert.m_p_MeV
    m_n=unitconvert.m_n_MeV
    n_s_mev4=unitconvert.toMev4(baryon_density_s,'mevfm')
    N1,N2=N1,N2
    crust_eos_array=crust_eos_array
    def chi(self,x):
        return (x*(1+x**2)**0.5*(2*x**2+1)-np.log(x+(1+x**2)**0.5))/(8*np.pi**2)
    def k_FQ(self,k_Fn):
        return (np.sign(9*k_Fn**2-self.kappa*self.Lambda*k_Fn-9*self.Lambda**2)+1)*((9*k_Fn**2-self.kappa*self.Lambda*k_Fn-9*self.Lambda**2)/(9*2*3*np.where(k_Fn>0,k_Fn,np.inf)))
    def k_FQ_jac(self,k_Fn):
        return (np.sign(9*k_Fn**2-self.kappa*self.Lambda*k_Fn-9*self.Lambda**2)+1)*((k_Fn**2+self.Lambda**2)/(2*3*np.where(k_Fn>0,k_Fn**2,np.inf)))
    def dndkF(self,k_Fn):
        return 3*k_Fn**2-81*self.k_FQ(k_Fn)**2*self.k_FQ_jac(k_Fn)
    def n_n(self,k_Fn):
        return (k_Fn**3-27*self.k_FQ(k_Fn)**3)/(3*np.pi**2)
    def Energy_density_n(self,k_Fn):
        k_0n=3*self.k_FQ(k_Fn)
        energy_density_kin=self.m_n**4*(self.chi(k_Fn/self.m_n)-self.chi(k_0n/self.m_n))
        n_n_mev4=self.n_n(k_Fn)
        n_n=unitconvert.toMevfm(n_n_mev4,'mev4')
        energy_density_pot=unitconvert.toMev4(self.Potential.eosDensity_from_n(n_n),'mevfm')
        return energy_density_kin+energy_density_pot

    def Chempo_PNM(self,k_Fn,n_n):
        mu_kin_n=(k_Fn**2+self.m_n**2)**0.5
        mu_pot_n=self.Potential.eosChempo_from_n(n_n)
        return mu_kin_n+mu_pot_n

    def Chempo(self,k_Fn):
        k_0n=3*self.k_FQ(k_Fn)
        Kn=(k_0n/k_Fn)**2*(1+(self.Lambda/k_Fn)**2)
        ((k_Fn**2+self.m_n**2)**0.5-Kn*(k_0n**2+self.m_n**2)**0.5)/(1-Kn)
        mu_kin_n=((k_Fn**2+self.m_n**2)**0.5-Kn*(k_0n**2+self.m_n**2)**0.5)/(1-Kn)
        n_n_mev4=self.n_n(k_Fn)
        n_n=unitconvert.toMevfm(n_n_mev4,'mev4')
        mu_pot_n=self.Potential.eosChempo_from_n(n_n)
        return mu_kin_n+mu_pot_n
    
    def Energy_density(self,k_Fn,mu_n):
        k_Fu=self.k_Fu(mu_n)
        k_Fd=2**(1/3)*k_Fu
        return self.Energy_density_n(k_Fn)+3*self.m_u**4*self.chi(k_Fu/self.m_u)+3*self.m_d**4*self.chi(k_Fd/self.m_d)
    def Pressure(self,k_Fn,mu_n):
        return mu_n*self.n_B(k_Fn,mu_n)-self.Energy_density(k_Fn,mu_n)
    def k_Fu(self,mu_n):
        mu_square=mu_n**2
        m1_square,m2_square=self.m_u**2,self.m_d**2
        tmp_213=2**(1/3)
        tmp_223=2**(2/3)
        k_Fu=np.sqrt((tmp_223*(4*m1_square-16*m2_square+4*mu_square)-m1_square+4*m2_square+mu_square-4*np.sqrt(8*tmp_213*m1_square*mu_square+tmp_223*(-m1_square*mu_square-4*m2_square*mu_square+mu_square**2)+m2_square*mu_square))/(-8*tmp_223+32*tmp_213+1))
        return k_Fu

    def n_B(self,k_Fn,mu_n):
        k_Fu=self.k_Fu(mu_n)
        return self.n_n(k_Fn)+k_Fu**3/(np.pi**2)  #2k_Fu**3=k_Fd**3, and g=2*3 spin*color
    
    def n_B_for_newton(self,k_Fn_modify,n_B):
        k_Fn=self.k_Fn_max-np.exp(k_Fn_modify)
        mu_n=self.Chempo(k_Fn)
        return self.n_B(k_Fn,mu_n)-n_B
    
    def __init__(self,args,defaut_u_max=12,potential=False,s_k=[0,2]):
        self.Lambda,self.n_t=args
        self.Potential=potential
        self.args=args
        self.k_F0,self.k_Ft=(unitconvert.toMev4(np.array([self.crust_eos_array[0][-1],self.n_t]),'mevfm')*(3*np.pi**2))**(1/3)
        self.mu_tn=self.Chempo_PNM(self.k_Ft,self.n_t)
        self.m_u=self.mu_tn/3
        self.m_d=self.mu_tn/3
        self.mass_args=(self.m_p,self.m_n,self.m_u,self.m_d)
        self.kappa=9*(self.k_Ft**2-self.Lambda**2)/(self.k_Ft*self.Lambda)
        if(self.kappa<0):
            for init_i in [500,800,1.1*self.k_Ft,300]:
                sol=root(self.dndkF,init_i,method='hybr')
                if(sol.success):
                    self.k_Fn_max=np.where(sol.x[0]<2000,sol.x[0],2000)
                    self.init_kFn_max=init_i
                    break
                self.k_Fn_max=2000
        else:
            self.k_Fn_max=2000
        for init_i in [0.01,1,2]:
            sol=root(self.n_B_for_newton,init_i,args=(self.n_s_mev4*defaut_u_max,))
            if(sol.success):
                self.init_nB_max=init_i
                break
        if(sol.success):
            self.k_Fmax=self.k_Fn_max-np.exp(sol.x[0])
            k_Fn_array=np.linspace(self.k_Ft,self.k_Fmax,self.N2)[1:]
            k_0n_array=3*self.k_FQ(k_Fn_array)
            mu_n_array=self.Chempo(k_Fn_array)
            k_Fu_array=self.k_Fu(mu_n_array)
            n_array  =self.n_B(k_Fn_array,mu_n_array)
            e_array  =self.Energy_density(k_Fn_array,mu_n_array)
            p_array  =mu_n_array*n_array-e_array
            k_Fn_array_PNM=np.linspace(1.1*self.k_F0,self.k_Ft,self.N1) #1.1 is to have a gap between crustal EOS and PNM EOS
            n_array_PNM= k_Fn_array_PNM**3/(3*np.pi**2)
            e_array_PNM= self.m_n**4*self.chi(k_Fn_array_PNM/self.m_n)+unitconvert.toMev4(self.Potential.eosDensity_from_n(unitconvert.toMevfm(n_array_PNM,'mev4')),'mevfm')
            p_array_PNM= ((self.m_n**2+k_Fn_array_PNM**2)**0.5+self.Potential.eosChempo_from_n(unitconvert.toMevfm(n_array_PNM,'mev4')))*n_array_PNM-e_array_PNM
            self.k_Fn_array=np.concatenate((k_Fn_array_PNM,k_Fn_array))
            self.k_0n_array=np.concatenate((0*k_Fn_array_PNM,k_0n_array))
            self.k_Fu_array=np.concatenate((0*k_Fn_array_PNM,k_Fu_array))
            self.eos_array_high=unitconvert.toMevfm(np.concatenate((np.array([n_array_PNM,e_array_PNM,p_array_PNM]),np.array([n_array,e_array,p_array])),axis=1),'mev4')
            self.eos_array=np.concatenate((self.crust_eos_array,self.eos_array_high),axis=1)
            EOS_interpolation.__init__(self,0.16,self.eos_array,s_k=s_k)
        else:
            print('Error!',args)
            self.eos_array=0
