#This is the script to control the trajectory plot output 
# of the simulation code. No changes are made to any of
# the fields of the simulation. Created by Debdut Sengupta 
# on 23 July 2022.

function plotTrajectory(sol,p)


    #Settings
    themeColour=:dark
    saveFigure=false
    zoomSize="large"
    legendLocation=:topleft

    #Attributes
    AU=p[5]
    ephemerides=p[33]
    colours=p[35]
    labels=p[36]
    if(zoomSize=="large")
        xlims=(-40.0,40.0)
        ylims=(-40.0,40.0)
        zlims=(-5.0,30.0)
    elseif(zoomSize=="medium")
        xlims=(-5.5,5.5)
        ylims=(-5.5,5.5)
        zlims=(-5.0,30.0)
    elseif(zoomSize=="small")
        xlims=(-2.0,2.0)
        ylims=(-2.0,2.0)
        zlims=(-5.0,10.0)
    end

    #Plotting sail
    theme(themeColour)
    plot(sol[1,:]./AU,sol[3,:]./AU,sol[5,:]./AU, 
        vars=(1, 3, 5), linecolor=:gold3, label="Sail", 
        xlabel="x (AU)", ylabel="y (AU)", zlabel="z (AU)", 
        plot_title="Trajectory of sail", xlims=xlims, 
        ylims=ylims, zlims=zlims, legend=legendLocation)

    #Plotting more bodies
    for k in eachindex(ephemerides)
        global eph = ephemerides[k]

        if (k==length(ephemerides))
            display(plot!(eph[:, 2]./AU, eph[:, 3]./AU, 
                eph[:, 4]./AU, linecolor=colours[k], label=labels[k]))
        else
            plot!(eph[:, 2]./AU, eph[:, 3]./AU, 
                eph[:, 4]./AU, linecolor=colours[k], label=labels[k])
        end
    end

    if(saveFigure)
        savefig("TrajectoryPlot")
    end

end
