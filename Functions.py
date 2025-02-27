# functions for MacAyeal_model_extended.py are defined in this file

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import math as m
import pandas as pd


def ODE(xx, epsilon, t, al, Q):
    dxdt = xx*al - epsilon*(xx**3) + Q
    
    return dxdt
    
def RK4(x, epsilon, t, alpha, Q, dt):
    
    dx = ODE(x, epsilon, t, alpha, Q); x1 = x + dx*dt/2
    dx1 = ODE(x1, epsilon, t, alpha, Q); x2 = x + dx1*dt/2
    dx2 = ODE(x2, epsilon, t, alpha, Q); x3 = x + dt*dx2
    dx2 = dx2+dx1
    dx1 = ODE(x3, epsilon, t, alpha, Q)
    
    return x + dx*dt/6 + dx2*dt/3 + dx1*dt/6

def Alphafun(alpha_1, alpha_2, lamb, A_1, A_2, R, switch, t):
    if switch==1: # alpha increases
        decay_con = lambd_a; alpha_lim = alpha_1; A = A_1
    else: # alpha decreases
        decay_con = lambd_a*R; alpha_lim = alpha_2; A = A_2
    
    alpha = ((A*m.exp(-1*decay_con*t))+alpha_lim)
    
    return alpha

def Solarforcing(q1, q0, omega_q, t):
    q = q1*m.sin(omega_q*t) + q0
    
    return q

def SolarModulator(A_m, qmag, qshift, omega_q, t):
    mod = A_m*m.sin(2*m.pi*t/100)
    AM = (qmag+mod)*m.sin(omega_q*t) + qshift
    
    return AM
    
def combined_forcing(A_m, q1, q0, omega_ob, omega_pr, t, c1, c2):
    mod = A_m*m.sin(2*m.pi*t/100)
    AM = c2*(q1+mod)*m.sin(omega_pr*t)
    q_ob = c1*q1*m.sin(omega_ob*t) 
    q = q_ob + AM + q0
    
    return q


def CuspLocus(alpha, epsilon):
    Qcusp = m.sqrt(4*(alpha**3)/(27*epsilon))
    
    return Qcusp
    
def Ramp(f1, f2, ramptime, t_ramp):
    f_ramp = ((-1*(f1-f2))/ramptime)*t_ramp + f1
    
    return f_ramp



