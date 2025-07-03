# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 22:01:17 2021

@author: HP
edited by: GBR on Tue Sep 28 10:46 2021
"""

import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# constants all in SI units
G = 6.6743e-11 #Gravity constant
M = 1.98847e30 #Sun mass
c = 2.99792458e8 #Speed of light
S = 1357 #Solar constant
AU = 1.495978707e11 #Astronomical unit

#sail properties
R = 5 #half of side
A = math.pi*R**2
d = 0.9 #surface density
C_spec = [0.719, 0.719] #specular reflectivity coefficient
C_abs = [0.162, 0.162] #absorption coefficient
C_diff = [0.117, 0.117] #diffuse reflectivity coefficient
B_f = 0 #Lambertian coefficient
K = 0 #thermal emissivity coefficient
h = 0 #center of pressure displacemnt
n = 0.7*math.pi/180  #membrane displacement angle
ksi = 0.2*math.pi/180 #membrane torsion angle
m = A*d #mass of membrane
I_s = 1/2*m*R**2 #moment of inertia about spin axis

#constant terms in dynamic equations
u = G*M #constant term for gravity
C = math.pi*R**3/3*S*AU**2/c #constant term for solar pressure

C_a = S*A*AU**2/(c*m)
C_n1 = []
C_n2 = []
C_r = []

C_n1.insert(0, 2*C_spec[0]) #constant term for specular reflection
C_n2.insert(0, B_f * C_diff[0]+ K*C_abs[0]) #constant term for normal emission
C_r.insert(0, C_abs[0] + C_diff[0]) #constant term for absorption

C_n1.insert(1, 2*C_spec[1]) #constant term for specular reflection
C_n2.insert(1, B_f * C_diff[1] + K*C_abs[1]) #constant term for normal emission
C_r.insert(1, C_abs[1] + C_diff[1]) #constant term for absorption

#initial vector format [x,xdot,y,ydot,z,zdot,alpha,delta,omega]
x0 = [-0.71843*AU,0,0,35.26e3,0,0,-40*math.pi/180,50*math.pi/180,-4.7/60*math.pi]
time = 10/12*31556926 #time of simulation

#derivative of ODE
def dx_dt(t,x):
    r = np.array([x[0],x[2],x[4]]) #radius of orbit
    r_mag = np.linalg.norm(r) #magnitude of radius
    r_norm = r/r_mag
    alpha_s = math.atan2(-x[2], -x[0])
    delta_s = math.asin(-x[4]/r_mag)
    norm = np.array([math.cos(x[7])*math.cos(x[6]),math.cos(x[7])*math.sin(x[6]),math.sin(x[7])])
    cosbeta = math.cos(alpha_s-x[6])*math.cos(delta_s-x[7])
    if(cosbeta > 0): #seting reflective side
        a = (-u*r_norm\
            + C_a*(C_r[0]*abs(-cosbeta)*r_norm\
            - (2/3*C_diff[0]+K*C_abs[0]+2*C_spec[0]*abs(cosbeta))*cosbeta*norm))/r_mag**2 #acceleration
        v = [-(2*C_n1[0]*abs(-cosbeta)+C_n2[0]+C_r[0])*ksi,
             -(2*C_n1[0]*abs(-cosbeta)+C_n2[0])*n-C_r[0]*cosbeta*(n+3*h/R),
             (2*C_n1[0]*cosbeta*abs(-cosbeta)+2*C_n2[0]*cosbeta+C_r[0]*(1-cosbeta**2))*ksi]
    else:
        a = (-u*r_norm\
            + C_a*(C_r[1]*abs(-cosbeta)*r_norm\
            -(2/3*C_diff[1]+K*C_abs[1]+2*C_spec[1]*abs(-cosbeta))*cosbeta*norm))/r_mag**2 #acceleration
        v = [-(2*C_n1[1]*abs(-cosbeta)+C_n2[1]+C_r[1])*ksi,
             -(2*C_n1[1]*abs(-cosbeta)+C_n2[1])*n-C_r[1]*abs(-cosbeta)*(n+3*h/R),
             (2*C_n1[1]*cosbeta*abs(-cosbeta)+2*C_n2[1]*cosbeta+C_r[1]*(1-cosbeta**2))*ksi]
    A = np.array([[1/(math.cos(x[7])*I_s*x[8]),0,0],
                  [0,1/(I_s*x[8]),0],
                  [0,0,1/(I_s)]])
    B = np.array([[math.cos(delta_s-x[7])*math.sin(alpha_s-x[6]),-math.sin(delta_s-x[7]),0],
                  [math.sin(delta_s-x[7]),math.cos(delta_s-x[7])*math.sin(alpha_s-x[6]),0],
                  [0,0,1]])
    e = np.dot(np.dot(A,B),v)*C/r_mag**2
    print(t/time)
    return(x[1],a[0],x[3],a[1],x[5],a[2],e[0],e[1],e[2]) #return vector

#ODE setup
t_span = [0, time] #input time
rtol = 1e-7 #relative tolerance
atol = 1e-9 #absolute tolerance
method = 'BDF' #method of solving, if RK45 fails try BDF
sol = solve_ivp(dx_dt, t_span, x0, rtol=rtol, atol=atol, method=method) #solver

#saving
np.savez_compressed('result',t=sol.t,y=sol.y)
