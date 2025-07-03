#This is the script to control the energy calculations as well
# as the energy plot output. Currently takes into account the
# Kinetic Energy and the gravitational potential energy wrt
# all bodies. Relativistic effects are not included and the srp
# energy is not included. The benefit of having this as a 
# separate file is that it can easily be switched on and off, 
# although it takes slightly longer to run. Created by Debdut 
# Sengupta on 11 September 2022.


include("time.jl")

function calcEnergyVec(sol,p)

    #Attributes
    ephemerides = p[33]
    masses      = p[34]
    m           = p[20]

    #Energy calculations
    KE          = zeros(size(sol.u))
    PE          = zeros(size(sol.u))

    for i in eachindex(sol.u)
        state = sol.u[i]

        #Kinetic energy, mv^2/2 (assuming no relativistic effects - quick check)
        KE[i] = state[2]^2+state[4]^2+state[6]^2

        ephTime = tToEphTime(sol.t[i])

        r = sqrt(state[1]^2+state[3]^2+state[5]^2)

        #Potential energy, sum(-GmM/r)

        for k in eachindex(ephemerides)

            r_i = r.-findState(ephemerides[k], ephTime)[1:3];

            r_mag_i = sqrt(r_i[1]^2+r_i[2]^2+r_i[3]^2);

            PE[i]-= p[1]*masses[k]*m/r_mag_i;

        end

    end
    KE.=KE.*(m/(2)) #KE=v^2 * m/2

    #Total Energy
    TE = KE.+PE

    return [KE,PE,TE]
end

function calcEnergy(t, x,vel, p)

    #WIP

    #Attributes
    ephemerides = p[33]
    masses      = p[34]
    m           = p[20]

    #Energy calculations

    #Kinetic energy, mv^2/2 (assuming no relativistic effects - quick check)
    KE    = m/2*(vel[1]^2+vel[2]^2+vel[3]^2)
    p[37] = [p[37];KE]

    ephTime = tToEphTime(t)

    PE = 0
    r = sqrt(x[1]^2+x[2]^2+x[3]^2)

    #Potential energy, sum(-GmM/r)

    for k in eachindex(ephemerides)

        r_ij = r.-findState(ephemerides[k], ephTime)[1:3];

        r_mag_ij = sqrt(r_ij[1]^2+r_ij[2]^2+r_ij[3]^2);

        PE-= p[1]*masses[k]*m/r_mag_ij;
    end

    p[38] = [p[38];PE]

    #Total Energy
    TE    = KE+PE
    p[39] = [p[39];TE]

    return [KE,PE,TE]
end

function plotEnergy(sol,p)

    if useCallback
        energies = calcEnergyVec(sol,p)
    end

    KE = energies[1]
    PE = energies[2]
    TE = energies[3]

    #Settings
    themeColour = :default
    legendLocation = :bottomright
    saveFigure = false

    theme(themeColour)

    #Converting energy to GJ and plotting
    plot(sol.t./31556926.0,KE./(1e9),linecolor="red", label="Kinetic Energy", 
        xlabel="Time (yrs)", ylabel="Energy (GJ)", plot_title="Energy of the sail 
        throughout the mission", legend=legendLocation)
    plot!(sol.t./31556926.0,PE./(1e9),linecolor="blue", label="Potential Energy")
    display(plot!(sol.t./31556926.0,TE./(1e9),linecolor="green",label="Total Energy"))

    if(saveFigure)
        savefig("EnergyPlot")
    end

end
