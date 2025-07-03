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

function gravity!(a, p, t, r)

    M = p[2]
    r_mag = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    r_norm = r / r_mag
    a .-= p[1] * M * r_norm / r_mag^2

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


function drag!(a, p, vel, r_mag)

    #normalised velocity vector
    v_mag = norm(vel)
    v_norm = vel / v_mag
    
    #constants
    p0 = 10^5 # pressure on the surface of the Earth in Pa
    MolarMass = 0.0289652 # molar mass of dry air in kg/mol
    R = 8.31446 #ideal gas constant in J/(mol·K)
    L = 0.0065 #temperature lapse rate in K/m
    g = 9.80665 #earth-surface gravitational acceleration in m/s2
    h0 = 6371*10^3 #radius of earth
    Cd = 1.05 #drag coefficient of a cube
    T0 = 15.04 + 273.1 #temperature at the surface in K

    #altitude above sea level
    h = r_mag - h0
    
    #atmosphere model: https://www.grc.nasa.gov/www/k-12/rocket/atmosmet.html
    if h <= 11000 #troposphere
        #temperature
        T = T0 - L * h
        #pressure
        pressure = p0 * (T / T0)^(g * MolarMass / (L * R))

    elseif h <= 25000 && h > 11000 #lower stratosphere
        #temperature
        T = 273.1 - 56.46
        #pressure
        pressure = 22.65 * exp(1.73 - (0.000157*h))
    
    else #upper stratosphere
        #temperature
        T = 273.1 - 131.21 + (0.00299 * h)
        #pressure
        pressure = 2.488 * (T / 216.6)^(-11.388)
    end

    #density
    density = (pressure * MolarMass / (R * T))*500
    

    #Uncomment the following if you want to see it spiralling inwards
    #density = 1e-11

    #update acceleration
    a .-= Cd * (v_mag^2) * density * p[7] * v_norm / (2 * p[20])

    #Comments:
    #Fix 1: Added . before -= in line 148 above for vector notation
    #Fix 2: Added the call to the drag function in dx_dt!.jl right after the gravity function (since having it below the change allocations means it doesn't affect anything)
    #Fix 3: Changed p[6] to p[7] in line 148 above since we need the area.
    #Fix 4: Added T0 (you had just missed it)
end



