#=
functions documentation
=#

"""
# gravity!(a, p, t, r, ephemerides)
updates the acceleration of the body due to gravity
## current state
only accounts for the Sun
## future possiblities
add more planets
"""
function gravity! end

"""
# srp!(a,p,r_norm,r_mag,cosbeta,norm,x,I_s,delta_s,alpha_s,C)
updates the acceleration of the body due to SRP effects
## current state
assumes flat sail with no deformations
## future possiblities
include deformations. develop the geometrical model. add effects from the spacecraft hub
"""
function srp! end

"""
# plasmaDrag!(a)
updates the acceleration of the body due to plasma drag
## current state
no effect on acceleration
## future possiblities
use data and modelling from Voyager missions on plasma density in the solar system to predict the effect of plasma drag
"""
function plasmaDrag! end

"""
# relativity!(a,p,r,vel,R_BIF,r_mag,vel_mag)
applies relativistic corrections to the acceleration of the body due to the orbit around the Sun
## current state
accounts for the Schwarzschild solution and Lense-Thirring Precession (the two corrections accounted for in GMAT for orbits around the Sun)
## future possiblities
R_BIF needs to be updated in main.jl
add more corrections
"""
function relativity! end

"""
# update_v!(p,cosbeta,x,I_s,delta_s,alpha_s,C,r_mag)
updates the angular velocities of the sail from the solar radiation pressure (called from srp!)
## current state
updates 3D angular velocity vector from solar radiation pressure
## future possiblities
update angular velocity vector from other forces as well
"""
function update_v! end

#=
functions implementation
=#

include("time.jl")

function quat_mul(q, q_prim)
    A = [q_prim[4] q_prim[3] -q_prim[2] q_prim[1]
        -q_prim[3] q_prim[4] q_prim[1] q_prim[2]
        q_prim[2] -q_prim[1] q_prim[4] q_prim[3]
        -q_prim[1] -q_prim[2] -q_prim[3] q_prim[4]]
    return A * q
end

function quat_inv(q)
    return [-q[1], -q[2], -q[3], q[4]]
end

function gravity!(a, p, t, r)

    ephemerides   = p[33]
    masses        = p[34]
    ephTime       = tToEphTime(t)
    for i in eachindex(ephemerides)
        r_i      = r .- findState(ephemerides[i], ephTime)[1:3]
        r_mag_i  = sqrt(r_i[1]^2 + r_i[2]^2 + r_i[3]^2)
        r_norm_i = r_i / r_mag_i
        a      .-= p[1] * masses[i] * r_norm_i / r_mag_i^2
    end


end

function srp!(a, e, p, r_norm, r_mag, cosbeta, norm, x, I_s, delta_s, alpha_s, C)
    #The distance to the Sun is no longer r, needs to be updated
    #Also needs to be updated with a more accurate model for when the sail is close to the Sun (no longer point source)


    if (cosbeta > 0)
        a .+= (p[24] * (p[27] * abs(-cosbeta) * r_norm
                        -
                        (2 / 3 * p[13] + p[16] * p[11] + 2 * p[9] * abs(cosbeta)) * cosbeta * norm)) / r_mag^2
        pressure = (p[24] * (p[27] * abs(-cosbeta) * r_norm
                             -
                             (2 / 3 * p[13] + p[16] * p[11] + 2 * p[9] * abs(cosbeta)) * cosbeta * norm)) / r_mag^2 * (p[30] / p[7])
    else
        a .+= (p[24] * (p[30] * abs(-cosbeta) * r_norm
                        -
                        (2 / 3 * p[14] + p[16] * p[12] + 2 * p[10] * abs(-cosbeta)) * cosbeta * norm)) / r_mag^2
        pressure = (p[24] * (p[30] * abs(-cosbeta) * r_norm
                             -
                             (2 / 3 * p[14] + p[16] * p[12] + 2 * p[10] * abs(-cosbeta)) * cosbeta * norm)) / r_mag^2 * (p[30] / p[7])
    end

    #Attitude update
    e .= update_v!(p, cosbeta, x, I_s, delta_s, alpha_s, C, r_mag, pressure)
    #e.=update_v_quaternion!(p,cosbeta,x,I_s,delta_s,alpha_s,C,r_mag,pressure)

end

function plasmaDrag!(a)
    a .+= 0
end

function gravityWithGR!(a, p, t, r, vel)

    ephemerides = p[33]
    β           = 1
    γ           = 1
    c           = p[3]
    masses      = p[34]
    G           = p[1]

    ephTime     = tToEphTime(t)
    vel_mag     = sqrt(vel[1]^2 + vel[2]^2 + vel[3]^2)

    for j in eachindex(ephemerides)
        r_ij = r .- findState(ephemerides[j], ephTime)[1:3]
        r_ij_mag = sqrt(r_ij[1]^2 + r_ij[2]^2 + r_ij[3]^2)
        vel_j = findState(ephemerides[j], ephTime)[4:6]
        vel_j_mag = sqrt(vel_j[1]^2 + vel_j[2]^2 + vel_j[3]^2) #Working now
        a_j = [0.0,0.0,0.0] #currently: need numeric differentiation from velocities. Can embed directly into file.

        term12 = 0.0
        term13 = 0.0
        for k in eachindex(ephemerides)

            r_ik = r .- findState(ephemerides[k], ephTime)[1:3]
            r_ik_mag = sqrt(r_ik[1]^2 + r_ik[2]^2 + r_ik[3]^2)
            term12 += G*masses[k]/r_ik_mag

            if(k!=j)
                r_jk = findState(ephemerides[j], ephTime)[1:3] .- findState(ephemerides[k], ephTime)[1:3]
                r_jk_mag = sqrt(r_jk[1]^2 + r_jk[2]^2 + r_jk[3]^2)
                term13 += G*masses[k]/r_jk_mag
            end

        end

        term12 *= -2*(β + γ)/(c^2)
        term13 *= -(2*β - 1)/(c^2)
        term14  = γ*(vel_mag/c)^2
        term15  = (1 + γ)*(vel_j_mag/c)^2
        term16  = 2*(1 + γ)/(c^2)*dot(vel,vel_j)
        term17  = -3/(2*c^2)*(dot(r_ij,vel_j)/r_ij_mag)^2
        term18  = 1/(2*c^2)*dot(-r_ij,a_j)

        #term 1
        a .-= G * masses[j] * r_ij / r_ij_mag^3 * (1 +
            term12 + term13 + term14 + term15 + term16 + term17 + term18)

        #term 2
        a .+= 1/(c^2)*G*masses[j]/r_ij_mag^3*
            dot(r_ij,(2 + 2*γ)*vel - (1 + 2*γ)*vel_j) * (vel - vel_j)

        #term 3
        a .+= (3 + 4*γ)/(2*c^2)*G*masses[j]*a_j/r_ij_mag
        
        
        #No Lense-Thirring Precession? Probably not necessary.
    end


end


function relativityOLD!(a, p, r, t, vel)

    #Can use cross(a,b) for a cross b

    #Old

    ##Add relativistic corrections from other planets

    #J = R_BIF * [0; 0; 0.4 * p[31]^2 * p[32]]

    ##Schwarzschild solution
    #a .+= p[22] / (p[3]^2 * r_mag^3) * ((4 * p[22] / r_mag - vel_mag^2) * r + 4 * (r' * vel) * vel)

    ##Lense-Thirring Precession
    #a .+= p[22] / (p[3]^2 * r_mag^3) * (3 / r_mag^2 * cross(r, vel) * (r' * J) + cross(vel, J))

end

function update_v!(p, cosbeta, x, I_s, delta_s, alpha_s, C, r_mag, pressure)

    n = norm(pressure) / (1420.0 * 7.5e-6 * 21.37 / 20 * x[9]^2) * (3 - 0.34) / 2 * 1.31199
    ksi = 0.02 * pi / 180

    if (cosbeta > 0) #setting reflective side
        v = [-(2 * p[25] * abs(-cosbeta) + p[26] + p[27]) * ksi,
            -(2 * p[25] * abs(-cosbeta) + p[26]) * n - p[27] * cosbeta * (n + 3 * p[17] / p[6]),
            (2 * p[25] * cosbeta * abs(-cosbeta) + 2 * p[26] * cosbeta + p[27] * (1 - cosbeta^2)) * ksi]
    else
        v = [-(2 * p[28] * abs(-cosbeta) + p[29] + p[30]) * ksi,
            -(2 * p[28] * abs(-cosbeta) + p[29]) * (-n) - p[30] * abs(-cosbeta) * ((-n) + 3 * p[17] / p[6]),
            (2 * p[28] * cosbeta * abs(-cosbeta) + 2 * p[29] * cosbeta + p[30] * (1 - cosbeta^2)) * ksi]
    end

    A = [1/(cos(x[8])*I_s*x[9]) 0.0 0.0
        0.0 1/(I_s*x[9]) 0.0
        0.0 0.0 1/(I_s)]
    B = [cos(delta_s - x[8])*sin(alpha_s - x[7]) -sin(delta_s - x[8]) 0.0
        sin(delta_s - x[8]) cos(delta_s - x[8])*sin(alpha_s - x[7]) 0.0
        0.0 0.0 1.0]
    return A * B * v * p[23] / r_mag^2 #Changed C to p[23]

end

function update_v_quaternion!(p, cosbeta, x, I_s, delta_s, alpha_s, C, r_mag, pressure)

    n = norm(pressure) / (1420.0 * 7.5e-6 * 21.37 / 20 * x[11]^2) * (3 - 0.34) / 2 * 1.31199
    ksi = 0.02 * pi / 180

    if (cosbeta > 0) #setting reflective side
        v = [-(2 * p[25] * abs(-cosbeta) + p[26] + p[27]) * ksi,
            -(2 * p[25] * abs(-cosbeta) + p[26]) * n - p[27] * cosbeta * (n + 3 * p[17] / p[6]),
            (2 * p[25] * cosbeta * abs(-cosbeta) + 2 * p[26] * cosbeta + p[27] * (1 - cosbeta^2)) * ksi]
    else
        v = [-(2 * p[28] * abs(-cosbeta) + p[29] + p[30]) * ksi,
            -(2 * p[28] * abs(-cosbeta) + p[29]) * (-n) - p[30] * abs(-cosbeta) * ((-n) + 3 * p[17] / p[6]),
            (2 * p[28] * cosbeta * abs(-cosbeta) + 2 * p[29] * cosbeta + p[30] * (1 - cosbeta^2)) * ksi]
    end

    #A = [1/(cos(x[8])*I_s*x[9]) 0.0 0.0
    #     0.0 1/(I_s*x[9]) 0.0
    #     0.0 0.0 1/(I_s)]
    alpha = 2.0 * atan(x[8], x[7])
    delta = pi / 2.0 - 2.0 * atan(x[7], x[10])
    #println(delta)
    B = [cos(delta_s - delta)*sin(alpha_s - alpha) -sin(delta_s - delta) 0.0
        sin(delta_s - delta) cos(delta_s - delta)*sin(alpha_s - alpha) 0.0
        0.0 0.0 1.0]
    #return A*B*v*p[23]/r_mag^2 #Changed C to p[23]
    q = [x[7], x[8], x[9], x[10]]
    T = B * v * p[23] / r_mag^2
    w_x = T[2] / (I_s * x[11])
    w_y = -T[1] / (I_s * x[11])
    omega = T[3] / I_s
    U = [-w_x, -w_y, 0.0, 0.0]
    qdot = 0.5 * quat_mul(q, U)
    return [qdot[1], qdot[2], qdot[3], qdot[4], omega]
    
end
