#=
functions documentation
=#

"""
# dx_dt!(dx,x,p,t)
updates the body's acceleration, position, and Euler angles
##  current state
takes into account several sources of acceleration. more can be added when uncommenting lines on the function definition
## future possiblities
"""
function dx_dt! end

#=
functions implementation
=#

include("forces.jl")

function dx_dt!(dx, x, p, t)

    #Vector properties
    r = [x[1], x[3], x[5]] #position vector
    vel = [x[2], x[4], x[6]] #velocity vector
    r_mag = sqrt(r[1]^2 + r[2]^2 + r[3]^2) #magnitude of radius
    r_norm = r / r_mag
    alpha_s = atan(-x[3], -x[1])
    delta_s = asin(-x[5] / r_mag)
    norm = [cos(x[8]) * cos(x[7]), cos(x[8]) * sin(x[7]), sin(x[8])]
    #println(norm)
    cosbeta = cos(alpha_s - x[7]) * cos(delta_s - x[8])
    R_BIF = [0 0 0
        0 0 0
        0 0 0] #Needs to be populated

    #acceleration update
    a = [0.0, 0.0, 0.0]
    e = [0.0, 0.0, 0.0]
    gravity!(a, p, t, r)
    drag!(a, p, vel, r_mag)
    #srp!(a, e, p, r_norm, r_mag, cosbeta, norm, x, I_s, delta_s, alpha_s, C)

    dx[1] = x[2]
    dx[2] = a[1]
    dx[3] = x[4]
    dx[4] = a[2]
    dx[5] = x[6]
    dx[6] = a[3]
    dx[7] = e[1]
    dx[8] = e[2]
    dx[9] = e[3]

    #drag!(a, p, vel, r_mag)
end
