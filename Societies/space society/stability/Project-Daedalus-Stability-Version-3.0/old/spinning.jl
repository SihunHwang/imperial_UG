using DifferentialEquations

# constants all in SI units
G = 6.6743e-11 #Gravity constant
M = 1.98847e30 #Sun mass
c = 2.99792458e8 #Speed of light
S = 1357 #Solar constant
AU = 1.495978707e11 #Astronomical unit

#sail properties
R = 5.0 #half of side
A = pi*R^2
d = 0.9 #surface density
C_spec = [0.719, 0.719] #specular reflectivity coefficient
C_abs = [0.162, 0.162] #absorption coefficient
C_diff = [0.117, 0.117] #diffuse reflectivity coefficient
B_f = 0 #Lambertian coefficient
K = 0 #thermal emissivity coefficient
h = 0 #center of pressure displacemnt
n = 0.7*pi/180  #membrane displacement angle
ksi = 0.2*pi/180 #membrane torsion angle
m = A*d #mass of membrane
I_s = 1/2*m*R^2 #moment of inertia about spin axis

#constant terms in dynamic equations
u = G*M #constant term for gravity
C = pi*R^3/3*S*AU^2/c #constant term for solar pressure

C_a = S*A*AU^2/(c*m)

C_n1_r = 2*C_spec[1]
C_n2_r = B_f * C_diff[1]+ K*C_abs[1]
C_r = C_abs[1] + C_diff[1]

C_n1_b = 2*C_spec[2]
C_n2_b = B_f * C_diff[2]+ K*C_abs[2]
C_b = C_abs[2] + C_diff[2]

p = [G, M, c, S, AU, R, A, d, C_spec[1], C_spec[2], C_abs[1], C_abs[2], C_diff[1], C_diff[2], B_f, K, h, n, ksi, m, I_s, u, C, C_a, C_n1_r, C_n2_r, C_r, C_n1_b, C_n2_b, C_b]
#p= [1, 2, 3, 4, 5,  6, 7, 8, 9,         10,        11,       12        13,        14,        15,  16,17,18,19,  20,21,  22,23,24,  25,     26,     27,  28,     29,     30]
function dx_dt!(dx,x,p,t)

    r = [x[1],x[3],x[5]] #radius of orbit
    r_mag = sqrt(x[1]^2+x[3]^2+x[5]^2) #magnitude of radius
    r_norm = r/r_mag
    alpha_s = atan(-x[3], -x[1])
    delta_s = asin(-x[5]/r_mag)
    norm = [cos(x[8])*cos(x[7]),cos(x[8])*sin(x[7]),sin(x[8])]
    cosbeta = cos(alpha_s-x[7])*cos(delta_s-x[8])
    if (cosbeta > 0) #seting reflective side
        a = (-u*r_norm
            + p[24]*(p[27]*abs(-cosbeta)*r_norm
            - (2/3*p[13]+p[16]*p[11]+2*p[9]*abs(cosbeta))*cosbeta*norm))/r_mag^2 #acceleration
        v = [-(2*p[25]*abs(-cosbeta)+p[26]+p[27])*p[19],
             -(2*p[25]*abs(-cosbeta)+p[26])*p[18]-p[27]*cosbeta*(p[18]+3*p[17]/p[6]),
             (2*p[25]*cosbeta*abs(-cosbeta)+2*p[26]*cosbeta+p[27]*(1-cosbeta^2))*p[19]]
    else
        a = (-u*r_norm
            + p[24]*(p[30]*abs(-cosbeta)*r_norm
            -(2/3*p[14]+p[16]*p[12]+2*p[10]*abs(-cosbeta))*cosbeta*norm))/r_mag^2 #acceleration
        v = [-(2*p[28]*abs(-cosbeta)+p[29]+p[30])*ksi,
             -(2*p[28]*abs(-cosbeta)+p[29])*p[18]-p[30]*abs(-cosbeta)*(p[18]+3*p[17]/p[6]),
             (2*p[28]*cosbeta*abs(-cosbeta)+2*p[29]*cosbeta+p[30]*(1-cosbeta^2))*p[19]]
    end
    A = [1/(cos(x[8])*I_s*x[9]) 0.0 0.0
         0.0 1/(I_s*x[9]) 0.0
         0.0 0.0 1/(I_s)]
    B = [cos(delta_s-x[8])*sin(alpha_s-x[7]) -sin(delta_s-x[8]) 0.0
        sin(delta_s-x[8]) cos(delta_s-x[8])*sin(alpha_s-x[7]) 0.0
        0.0 0.0 1.0]
    e = A*B*v*C/r_mag^2
    #println(t)

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


x0 = [-0.71843*AU,0.0,0.0,35.26e3,0.0,0.0,-40*pi/180,50*pi/180,-4.7/60*pi]
tspan = (0.0,1.0*31556926.0)

prob = ODEProblem(dx_dt!,x0,tspan,p)
#RadauIIA5()
sol = solve(prob, FBDF(), reltol=1e-7, abstol=1e-9)

using Plots
plot(sol,vars=(1,3,5))