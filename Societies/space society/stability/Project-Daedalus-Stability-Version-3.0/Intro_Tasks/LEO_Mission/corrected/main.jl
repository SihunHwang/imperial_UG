#=
Written by DS
Work started on Mar 24th 2022
Developed to analyse the 6 DoF stability of the Daedalus solar sailcraft in LEO Mission
=#

# packages
using DifferentialEquations
using LinearAlgebra
using Plots
using ForwardDiff

# other scripts required
include("dx_dt!.jl")

# constants
G = 6.6743e-11 #gravitational constant [N m^2 kg^-2]
c = 2.99792458e8 #speed of light [m s^-1]
S = 1357 #solar constant [W m^-2]
AU = 1.495978707e11 #astronomical unit [km]

# central object (Earth) properties [kg]
M = 5.972e24 #Mass of the central body

#sail properties
R = 21.37 / 2.0 #half of side
A = pi * R^2 #area of sail
d = 0.009 #surface density (mass per unit area)
C_spec = [0.719, 0.545]#[0.719, 0.719] #specular reflectivity coefficient
C_abs = [0.162, 0.162]#[0.162, 0.162] #absorption coefficient
C_diff = [0.117, 0.466] #diffuse reflectivity coefficient
B_f = 0 #Lambertian coefficient
K = 0 #thermal emissivity coefficient
h = 0 #center of pressure displacemnt
n = 0.007 * pi / 180  #membrane displacement angle
ksi = 0.02 * pi / 180 #membrane torsion angle
m = A * d #mass of membrane
I_s = 1 / 2 * m * R^2 #moment of inertia about spin axis

#constant terms in dynamic equations
u = G * M #constant term for gravity
C = pi * R^3 / 3 * S * AU^2 / c #constant term for solar pressure

C_a = S * A * AU^2 / (c * m)

C_n1_r = 2 * C_spec[1]
C_n2_r = B_f * C_diff[1] + K * C_abs[1]
C_r = C_abs[1] + C_diff[1]

C_n1_b = 2 * C_spec[2]
C_n2_b = B_f * C_diff[2] + K * C_abs[2]
C_b = C_abs[2] + C_diff[2]

#Solver parameters
span_y = 0.05 #mission span in y
tspan = (0.0, span_y * 31556926.0)

p = [G, M, c, S, AU, R, A, d, C_spec[1], C_spec[2], C_abs[1], C_abs[2], C_diff[1], C_diff[2], B_f, K, h, n, ksi, m, I_s, u, C, C_a, C_n1_r, C_n2_r, C_r, C_n1_b, C_n2_b, C_b]

#p= [1, 2, 3, 4, 5,  6, 7, 8, 9,         10,        11,       12        13,        14,        15,  16,17,18,19,  20,21,  22,23,24,  25,     26,     27,  28,     29,     30,   31,    32  ,      33    ,   34  ,    35  ,   36  ]

x0      = [6871e3, 0.0, 0.0, 7.6e3, 0.0, 0.0, 30.0 * pi / 180, 0.0, 1.0] #Relative to Earth
#x0     = [x     , vx , y  , vy    , z , vz , R.A.           , DEC, spin rate]
prob = ODEProblem(dx_dt!, x0, tspan, p)

maxiters = Int(1e7)
#If the solver does not work when you add Drag, use DP5() instead
sol = solve(prob, DP5(), reltol=1e-15, abstol=1e-15, maxiters=maxiters)

theme(:dark)
plot(sol, idxs=(1,3,5))
