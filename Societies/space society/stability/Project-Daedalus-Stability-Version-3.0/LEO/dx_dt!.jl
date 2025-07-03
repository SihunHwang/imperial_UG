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

"""
# dx_dt_quaternion!(dx,x,p,t)
updates the body's acceleration, position, and Euler angles, using Quaternions
##  current state
takes into account several sources of acceleration. more can be added when uncommenting lines on the function definition. uses Quaternions rather than orthogonal coordinate approach
## future possiblities
 
"""
function dx_dt_quaternion! end

#=
functions implementation
=#

include("forces.jl")

function dx_dt!(dx, x, p, t)

    #Vector properties
    r       = [x[1], x[3], x[5]] #position vector
    vel     = [x[2], x[4], x[6]] #velocity vector
    r_mag   = sqrt(r[1]^2 + r[2]^2 + r[3]^2) #magnitude of radius
    r_norm  = r / r_mag
    alpha_s = atan(-x[3], -x[1])
    delta_s = asin(-x[5] / r_mag)
    norm    = [cos(x[8]) * cos(x[7]), cos(x[8]) * sin(x[7]), sin(x[8])]
    #println(norm)
    cosbeta = cos(alpha_s - x[7]) * cos(delta_s - x[8])
    #R_BIF   = [0 0 0
    #    0 0 0
    #    0 0 0] #Needs to be populated

    #acceleration update
    a = [0.0, 0.0, 0.0]
    e = [0.0, 0.0, 0.0]

    if includeGR
        gravityWithGR!(a, p, t, r, vel)
    else
        gravity!(a, p, t, r)
    end

    srp!(a, e, p, r_norm, r_mag, cosbeta, norm, x, I_s, delta_s, alpha_s, C)
    #plasmaDrag!(a)

    dx[1] = x[2]
    dx[2] = a[1]
    dx[3] = x[4]
    dx[4] = a[2]
    dx[5] = x[6]
    dx[6] = a[3]
    dx[7] = e[1]
    dx[8] = e[2]
    dx[9] = e[3]

end

function dx_dt_quaternion!(dx, x, p, t)

    #Vector properties
    r = [x[1], x[3], x[5]] #position vector
    vel = [x[2], x[4], x[6]] #velocity vector
    r_mag = sqrt(r[1]^2 + r[2]^2 + r[3]^2) #magnitude of radius
    vel_mag = sqrt(x[2]^2 + x[4]^2 + x[6]^2)
    r_norm = r / r_mag
    alpha_s = atan(-x[3], -x[1])
    delta_s = asin(-x[5] / r_mag)
    alpha = 2.0 * atan(x[8], x[7])
    delta = pi / 2.0 - 2.0 * atan(x[7], x[10])
    norm = [cos(delta) * cos(alpha), cos(delta) * sin(alpha), sin(delta)]
    #println(norm)
    cosbeta = cos(alpha_s - alpha) * cos(delta_s - delta)
    R_BIF = [0 0 0
        0 0 0
        0 0 0] #Needs to be populated


    #acceleration update
    a = [0.0, 0.0, 0.0]
    e = [0.0, 0.0, 0.0, 0.0, 0.0]
    gravity!(a, p, t, r, ephemerides)

    srp!(a, e, p, r_norm, r_mag, cosbeta, norm, x, I_s, delta_s, alpha_s, C)
    #plasmaDrag!(a)
    #relativity!(a,p,r,vel,R_BIF,r_mag,vel_mag)

    dx[1] = x[2]
    dx[2] = a[1]
    dx[3] = x[4]
    dx[4] = a[2]
    dx[5] = x[6]
    dx[6] = a[3]
    dx[7] = e[1]
    dx[8] = e[2]
    dx[9] = e[3]
    dx[10] = e[4]
    dx[11] = e[5]

end
