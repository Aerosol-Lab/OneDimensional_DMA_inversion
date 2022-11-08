import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

k_B = 1.380650524e-23 # J/K
q_e = 1.602e-19       # C

# Air density
def Rho_gas(T_g):
    P_g = 101300
    R = 287 #J/kg·K
    rho = P_g/(R * T_g)
    return rho

# Gas dinamic viscosity
def Mu_gas(T_g):
    mu_g = (18.203E-6)*(293.15+110)/(T_g+110)*(T_g/293.15)**(1.5)
    return mu_g

# Gas mean free path
def Lambda_gas(T_g):
    P_g = 101300
    lambda_g = 66.5E-9*(101300/P_g)*(T_g/293.15)*(1+110/293.15)/(1+110/T_g)
    return lambda_g

# Gas Knudsen number
def GET_Cc(Dp,T_g):
    A1 = 1.142
    A2 = 0.558
    A3 = 0.999
    lambda_g = Lambda_gas(T_g)
    Kn = lambda_g/(0.5*Dp)
    return 1+A1*Kn+A2*Kn*np.exp(-A3/Kn)

# Friction coefficient
def friction(Dp,T_g):
    mu_g = Mu_gas(T_g)
    Cc = GET_Cc(Dp,T_g)
    f =  3*np.pi*mu_g*Dp/Cc
    return f
