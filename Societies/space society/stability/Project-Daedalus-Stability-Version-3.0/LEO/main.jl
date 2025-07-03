#=
Written by GBR and DS
Work started on Mar 24th 2022
Based on work by Piotr over summer 2021 on Python in the \old folder
Developed to analyse the 6 DoF stability of the Daedalus solar sailcraft
=#

# To do: 
# 1. Add relativity. 
# 2. Make animations 
# 3. Make quaternions work
# 4. Validate against GMAT
# 5. Add plasma drag
# 6. Find optimal initial conditions
# 7. Investigate stability in different conditions
# 8. Implement Callback function

# packages
using DifferentialEquations
using LinearAlgebra
using Plots
using ForwardDiff

# other scripts required
include("dx_dt!.jl")
include("readEphemeris.jl")
include("plotTrajectory.jl")
include("energy.jl")
include("time.jl")

# constants
G             = 6.6743e-11 #gravitational constant [N m^2 kg^-2]
c             = 2.99792458e8 #speed of light [m s^-1]
S             = 1357 #solar constant [W m^-2]
AU            = 1.495978707e11 #astronomical unit [km]

#Positions and masses of all bodies, ordered by standard radial position
ephemerides   = loadEphemerides();
masses        = [1.9885e30, 3.302e23, 4.8685e24, 5.972e24, 6.4171e23, 1.898e27, 5.6834e26, 8.6813e25, 1.02409e26]; 
colours       = [:yellow, :grey83, :orange, :blue, :red, :tan4, :yellow3, :lightskyblue1, :dodgerblue2]
labels        = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

# central object (Sun) properties [kg]
M             = masses[1] #Mass of the central body

#Terms for the Lense-Thirring Precession, to be removed.
#R_B     = 695508000 #Radius of central body
#omega_B = 2 * pi / (27 * 24 * 3600) #spin rate, varies along longitude

#sail properties
R             = 21.37 / 2.0 #half of side
A             = pi * R^2 #area of sail
d             = 0.009 #surface density (mass per unit area)
C_spec        = [0.719, 0.545]#[0.719, 0.719] #specular reflectivity coefficient
C_abs         = [0.162, 0.162]#[0.162, 0.162] #absorption coefficient
C_diff        = [0.117, 0.466] #diffuse reflectivity coefficient
B_f           = 0 #Lambertian coefficient
K             = 0 #thermal emissivity coefficient
h             = 0 #center of pressure displacemnt
n             = 0.007 * pi / 180  #membrane displacement angle
ksi           = 0.02 * pi / 180 #membrane torsion angle
m             = A * d #mass of membrane
I_s           = 1 / 2 * m * R^2 #moment of inertia about spin axis

#constant terms in dynamic equations
u             = G * M #constant term for gravity
C             = pi * R^3 / 3 * S * AU^2 / c #constant term for solar pressure

C_a           = S * A * AU^2 / (c * m)

C_n1_r        = 2 * C_spec[1]
C_n2_r        = B_f * C_diff[1] + K * C_abs[1]
C_r           = C_abs[1] + C_diff[1]

C_n1_b        = 2 * C_spec[2]
C_n2_b        = B_f * C_diff[2] + K * C_abs[2]
C_b           = C_abs[2] + C_diff[2]

#Energy initialisations
KE            = [0.0]
PE            = [0.0]
TE            = [0.0]

#Solver parameters
span_y        = 40 #mission span in y
tspan         = (0.0, yrsToSeconds(span_y))
useQuaternion = false
useCallback   = true
includeGR     = false
maxiters      = Int(1e7)


p = [G, M, c, S, AU, R, A, d, C_spec[1], C_spec[2], C_abs[1], C_abs[2], C_diff[1], C_diff[2], B_f, K, h, n, ksi, m, I_s, u, C, C_a, C_n1_r, C_n2_r, C_r, C_n1_b, C_n2_b, C_b, 0,   0   , ephemerides, masses, colours, labels, KE, PE, TE]
#p= [1, 2, 3, 4, 5,  6, 7, 8, 9,         10,        11,       12        13,        14,        15,  16,17,18,19,  20,21,  22,23,24,  25,     26,     27,  28,     29,     30,  31,  32  ,    33    ,      34  ,   35  ,    36    , 37, 38, 39]

if (useQuaternion)
    x0_quat = [-AU, 0.0, 0.0, 29.78e3, 0.0, 0.0, 
              sin((pi / 2.0 - x0[8]) / 2.0) * cos(x0[7] / 2.0), 
              sin((pi / 2.0 - x0[8]) / 2.0) * sin(x0[7] / 2.0), 
              cos((pi / 2.0 - x0[8]) / 2.0) * sin(x0[7] / 2.0), 
              cos((pi / 2.0 - x0[8]) / 2.0) * cos(x0[7] / 2.0), 
              1.0]
    prob = ODEProblem(dx_dt_quaternion!, x0_quat, tspan, p)
else
    x0      = [-AU, 0.0, 0.0, 29.78e3, 0.0, 0.0, 30.0 * pi / 180, 0.0, 1.0] # -40*pi/180, 50*pi/180, -4.7/60*pi
    #x0     = [x  , vx , y  , vy     , z  , vz , R.A.           , DEC, spin rate]
    prob = ODEProblem(dx_dt!, x0, tspan, p)
end

if (useCallback)
    #Callback function
    function condition(u,t,integrator)

        x        = [u[1],u[3],u[5]]
        vel      = [u[2],u[4],u[6]]
        energies = calcEnergy(t,x,vel,p)
        PE       = energies[2]
        TE       = energies[3]

        return TE>0 && abs(PE)<2e8
    end
    affect!(integrator) = integrator.t = tspan[2]
    #affect!(integrator) = terminate!(integrator)
    cb                  = DiscreteCallback(condition, affect!)
    sol                 = solve(prob, DP5(), reltol=1e-15, abstol=1e-15, maxiters=maxiters, callback=cb)
else
    sol = solve(prob, DP5(), reltol=1e-15, abstol=1e-15, maxiters=maxiters)
    #Other good solvers:
    #RadauIIA5()
    #RadauIIA5(autodiff=false)
end

plotEnergy(sol,p)
plotTrajectory(sol,p)



#Plotting ideas-------------------------------------------------------------------------
#animate(sol,vars=(1,3,5),lw=3,every=10) #Animates the orbit in time
#display(plot(sol, vars = (0,7)))
#display(plot(sol, vars = (0,8)))
#display(plot(sol, vars = (0,9)))
#animate(sol,vars=(1,3,5),lw=3,every=10) #Animates the orbit in time


#BVP solver-----------------------------------------------------------------------------

#function bc!(residual, u, p, t)
#    residual[1] = u[1][1] + AU
#    residual[2] = u[1][3]
#    residual[3] = u[1][5]
#    residual[1] = u[end][1] + AU
#    residual[2] = u[end][3]
#    residual[3] = u[end][5]
#end

#u0 = x0
#opt = BVProblem(dx_dt!, bc!, u0, tspan, p)
#bvp_sol = solve(opt, Shooting(DP5()),force_dtmin=true,abstol=1e-13,reltol=1e-13)
#bvp_sol = solve(opt, MIRK4(), dt=tspan[2]/200.0)
#display(plot(bvp_sol, vars = (1,3,5)))
#display(plot(sol, vars = (0,1)))
#display(plot(sol, vars = (0,3)))
#display(plot(sol, vars = (0,5)))
#display(plot(sol, vars = (0,2)))
#display(plot(sol, vars = (0,4)))
#display(plot(sol, vars = (0,6)))
#display(plot(sol, vars = (0,7)))
#display(plot(sol, vars = (0,8)))
#display(plot(sol, vars = (0,9)))
