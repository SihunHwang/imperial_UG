# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 11:44:04 2022

@author: HP
"""

using DifferentialEquations
using Plots

function membrane_fun(dy, y, p, r)
    h = 7.5e-6
    E = 3.19e9
    v = 0.34
    rho = 1420
    alpha = 9.8
    omega = 1#15*2*math.pi
    angle = 90 * pi / 180
    sin_phi = y[1] / sqrt(y[1]^2 + y[2]^2)
    cos_phi = y[2] / sqrt(y[1]^2 + y[2]^2)
    sin_phi = ifelse(isnan(sin_phi), 0.0, sin_phi)#np.nan_to_num(sin_phi)
    cos_phi = ifelse(isnan(cos_phi), 0.0, cos_phi)#np.nan_to_num(cos_phi)
    #q_z = -0.000217#-rho*h*alpha
    q_r = rho * h * (r + y[4]) * omega^2
    q_n = 0.000217 * (sin(angle - asin(sin_phi)))^2
    q_z = -q_n * cos_phi
    #q_phi = q_r*cos_phi
    q_phi = q_r * cos_phi + q_z * sin_phi
    dOmega_dr = -r * q_z
    dPsi_dr = (y[3] - q_r * r^2) / r
    dGamma_dr = (y[2] + E * h * r * (cos_phi - 1) - v * q_phi * r^2) / r
    du_dr = cos_phi * (1 + (sqrt(y[1]^2 + y[2]^2) - v * y[3]) / (E * h * r)) - 1
    dw_dr = sin_phi * (1 + (sqrt(y[1]^2 + y[2]^2) - v * y[3]) / (E * h * r))

    dy[1] = dOmega_dr
    dy[2] = dPsi_dr
    dy[3] = dGamma_dr
    dy[4] = du_dr
    dy[5] = dw_dr
end

function membrane_fun_norm(dy, y, p, r)
    #h = 7.5e-6
    v = 0.34
    #rho = 1420
    #alpha = 9.8
    #b = 21.37/2.0
    #omega = 1.0#15*2*math.pi
    angle = 90 * pi / 180
    E = p[1]#3.19e9/(rho*b^2*omega^2)
    #q_n0 = p[2]
    sin_phi = y[1] / sqrt(y[1]^2 + y[2]^2)
    cos_phi = y[2] / sqrt(y[1]^2 + y[2]^2)
    sin_phi = ifelse(isnan(sin_phi), 0.0, sin_phi)#np.nan_to_num(sin_phi)
    cos_phi = ifelse(isnan(cos_phi), 0.0, cos_phi)#np.nan_to_num(cos_phi)
    #q_z = -0.000217#-rho*h*alpha
    q_r = r + y[4]#rho*h*(r+y[4])*omega^2
    q_n = p[2] * (sin(angle - asin(sin_phi)))^2#0.000217/(rho*h*b*omega^2)
    q_z = -q_n * cos_phi
    #q_phi = q_r*cos_phi
    q_phi = q_r * cos_phi + q_z * sin_phi
    dOmega_dr = -r * q_z
    dPsi_dr = (y[3] - q_r * r^2) / r
    dGamma_dr = (y[2] + E * r * (cos_phi - 1) - v * q_phi * r^2) / r
    du_dr = cos_phi * (1 + (sqrt(y[1]^2 + y[2]^2) - v * y[3]) / (E * r)) - 1
    dw_dr = sin_phi * (1 + (sqrt(y[1]^2 + y[2]^2) - v * y[3]) / (E * r))

    dy[1] = dOmega_dr
    dy[2] = dPsi_dr
    dy[3] = dGamma_dr
    dy[4] = du_dr
    dy[5] = dw_dr
end

function membrane_bc(residual, y, p, r)
    v = 0.34
    residual[1] = y[end][1]
    residual[2] = y[end][2]
    #residual[3] = v*sqrt(y[1][1]^2 + y[1][2]^2) - y[1][3]
    residual[4] = y[1][4]
    residual[5] = y[1][5]
    #println(y)
    #return np.array([yb[0], yb[1], v*np.sqrt(np.square(ya[0])+np.square(ya[1]))-ya[2], ya[3], ya[4]])
end

#=
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
=#

rho = 1420.0
b = 21.37 / 2.0
omega = 1.0
h = 7.5e-6
v = 0.34

tspan = (0.33 / b, 1.0)
u0 = [0.0, 10.0 * tspan[1], v * 10.0 * tspan[1], 0.0, 0.0]


for j in 1:20
    for i in 0:20
        global p = [j / 5.0 * 3.19e9 / (rho * b^2 * omega^2), i / 5.0 * 0.000217 / (rho * h * b * omega^2)]
        opt = BVProblem(membrane_fun_norm, membrane_bc, u0, tspan, p)
        #bvp_sol = solve(opt, Shooting(DP5()),force_dtmin=true,abstol=1e-9,reltol=1e-7)
        bvp_sol = solve(opt, GeneralMIRK4(), dt=tspan[2] / 20.0)
        #display(plot(bvp_sol, vars = (0,1)))
        #display(plot(bvp_sol, vars = (0,2)))
        display(plot(bvp_sol, vars=(0, 3)))
        display(plot(bvp_sol, vars=(0, 5)))
        if (bvp_sol[3, end] < 0.0)
            break
        end
        println(bvp_sol[5, end])
    end
    println("")
end
