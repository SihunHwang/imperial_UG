# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 11:44:04 2022

@author: HP
"""

import numpy as np
import math
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

def fun(x, y):
    return np.vstack((y[1], -np.exp(y[0])))

def bc(ya, yb):
    return np.array([ya[0], yb[0]])

def membrane_fun(r,y):
    h = 7.5e-6
    E = 3.19e9
    v = 0.34
    rho = 1420
    alpha = 9.8
    omega = 1#15*2*math.pi
    angle = 90*math.pi/180
    sin_phi = y[0]/np.sqrt(np.square(y[0])+np.square(y[1]))
    cos_phi = y[1]/np.sqrt(np.square(y[0])+np.square(y[1]))
    sin_phi = np.nan_to_num(sin_phi)
    cos_phi = np.nan_to_num(cos_phi)
    #q_z = -0.000217#-rho*h*alpha
    q_r = rho*h*(r+y[3])*omega**2
    q_n = 0.000217*np.square(np.sin(angle-np.arcsin(sin_phi)))
    q_z = -q_n*cos_phi
    q_phi = q_r*cos_phi
    #q_phi = q_r*cos_phi + q_z*sin_phi
    dOmega_dr = -r*q_z
    dPsi_dr = (y[2]-q_r*r**2)/r
    dGamma_dr = (y[1]+E*h*r*(cos_phi-1)-v*q_phi*r**2)/r
    du_dr = cos_phi*(1+(np.sqrt(np.square(y[0])+np.square(y[1]))-v*y[2])/(E*h*r))-1
    dw_dr = sin_phi*(1+(np.sqrt(np.square(y[0])+np.square(y[1]))-v*y[2])/(E*h*r))
    return np.vstack((dOmega_dr, dPsi_dr, dGamma_dr, du_dr, dw_dr))

def membrane_bc(ya, yb):
    v = 0.34
    return np.array([yb[0], yb[1], v*np.sqrt(np.square(ya[0])+np.square(ya[1]))-ya[2], ya[3], ya[4]])


r_a = 0.33
r_b = 21.37/2

r = np.linspace(r_a, r_b, 2)
y = np.zeros((5, r.size))

res = solve_bvp(membrane_fun, membrane_bc, r, y, tol=1e-1)

x = np.linspace(r_a, r_b, 100)
y = res.sol(x)

angle = np.arcsin(y[0]/np.sqrt(np.square(y[0])+np.square(y[1])))

plt.plot(x,y[4])
#plt.plot(x,angle*180/math.pi)
#plt.plot(x,y[2]/x/7.5e-6)
plt.xlabel("radius [m]")
plt.ylabel("deflection [m]")
plt.show()