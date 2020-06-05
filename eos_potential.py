#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 15:11:14 2020

@author: sotzee
"""

import numpy as np
from sympy import symbols, diff,lambdify
import unitconvert
from config import saturation_density

class Fermions(object):
    ns=saturation_density
    def __init__(self,args):
        self.name, self.m, self.g=args #m in unit MeV, g is degenracy in spin or isospin...
    def set_mass(self,mass):
        self.m=mass
    def chi(self,x):
        return self.g*(x*(1+x**2)**0.5*(2*x**2+1)-np.log(x+(1+x**2)**0.5))/(16*np.pi**2)
    def phi(self,x): #x=kF/m demensionless
        return self.g*(x*(1+x**2)**0.5*(2*x**2-3)+3*np.log(x+(1+x**2)**0.5))/(48*np.pi**2)
    def psi(self,x):
        return self.g*(4*x**5/(1+x**2)**0.5-3*x*(1+x**2)**0.5*(2*x**2-3)-9*np.log(x+(1+x**2)**0.5))/(72*np.pi**2)
    def eosDensity_from_x(self,x,x0=0):
        return unitconvert.toMevfm(self.m**4*(self.chi(x)-self.chi(x0)),'mev4')
    def eosPressure_from_x(self,x,x0=0):
        return unitconvert.toMevfm(self.m**4*(self.phi(x)-self.phi(x0)),'mev4')
    def eosN3d2Edn2_from_x(self,x,x0=0):
        return unitconvert.toMevfm(self.m**4*(self.psi(x)-self.psi(x0)),'mev4')
    def eosCs2(self,x):
        return (2*self.eosPressure_from_x(x)+self.eosN3d2Edn2_from_x(x))/(self.eosDensity_from_x(x)+self.eosPressure_from_x(x))
    def eosBaryonDensity_from_x(self,x,x0=0):
        return unitconvert.toMevfm(self.g*((x*self.m)**3-(x0*self.m)**3)/(6*np.pi**2),'mev4')
    def eosChempo_from_x(self,x):
        return self.m*(x**2+1)**0.5
    def eosX_from_n(self,n):
        return np.sign(n)*np.abs(unitconvert.toMev4(n,'mevfm')*(6*np.pi**2/(self.g*self.m**3)))**(1/3)

class Potential_single(object):
    ns=saturation_density
    def __init__(self,args,sym_list,mean_potential_expr):
        self.args=args
        args_sym_list=sym_list[:-1]
        mean_potential_expr_subs=mean_potential_expr.subs(zip(args_sym_list,args))
        self.mean_potential_E=lambdify(sym_list[-1],mean_potential_expr_subs)
        self.mean_potential_dEdn=lambdify(sym_list[-1],diff(mean_potential_expr_subs,sym_list[-1]))
        self.mean_potential_d2Edn2=lambdify(sym_list[-1],diff(mean_potential_expr_subs,sym_list[-1],2))
        self.E=self.mean_potential_E(self.ns)
        self.L=3*self.ns*self.mean_potential_dEdn(self.ns)
        self.K=9*self.ns**2*self.mean_potential_d2Edn2(self.ns)
    def eosDensity_from_n(self,n):
        return n*self.mean_potential_E(n)
    def eosPressure_from_n(self,n):
        return n**2*self.mean_potential_dEdn(n)
    def eosChempo_from_n(self,n):
        return (self.eosDensity_from_n(n)+self.eosPressure_from_n(n))/n

crust_core_density=0.4*saturation_density #Not used in quarkyonic EOS
proton=Fermions(['proton',unitconvert.m_p_MeV,2])
neutron=Fermions(['neutron',unitconvert.m_n_MeV,2])
n0ns=np.array([0.4*saturation_density,saturation_density])
xs_p_sym=proton.eosX_from_n(n0ns/2)
xs_n_sym=neutron.eosX_from_n(n0ns/2)
xs_pnm=neutron.eosX_from_n(n0ns)
E_kin_sym=(proton.eosDensity_from_x(xs_p_sym)+neutron.eosDensity_from_x(xs_n_sym))/n0ns
L_kin_sym=3*(proton.eosPressure_from_x(xs_p_sym)+neutron.eosPressure_from_x(xs_n_sym))/n0ns
K_kin_sym=9*(proton.eosN3d2Edn2_from_x(xs_p_sym)+neutron.eosN3d2Edn2_from_x(xs_n_sym))/n0ns
ELK_kin_sym=np.array([E_kin_sym,L_kin_sym,K_kin_sym])
E_kin_pnm=neutron.eosDensity_from_x(xs_pnm)/n0ns
L_kin_pnm=3*neutron.eosPressure_from_x(xs_pnm)/n0ns
K_kin_pnm=9*neutron.eosN3d2Edn2_from_x(xs_pnm)/n0ns
ELK_kin_pnm=np.array([E_kin_pnm,L_kin_pnm,K_kin_pnm])

def V_Lattimer(n_s,a,b,gamma,n):
    return a*(n/n_s)+b*(n/n_s)**gamma
def fit_lattimer_pnm(para,ELgamma):
    Potential_Lattimer_pnm=Potential_single(np.concatenate((para,[ELgamma[2]])),syms_Lattimer,V_Lattimer_expr)
    EL_potential_pnm=np.array([Potential_Lattimer_pnm.E,Potential_Lattimer_pnm.L])
    return ELK_kin_pnm[:2,1]+EL_potential_pnm-np.array(ELgamma[:2])

sym_a, sym_b, sym_d, sym_gamma, sym_alpha, sym_beta, sym_n= symbols('a b d gamma alpha beta n', real=True)
syms_Lattimer=[sym_a, sym_b, sym_gamma, sym_n]
V_Lattimer_expr=V_Lattimer(saturation_density, sym_a, sym_b, sym_gamma, sym_n)