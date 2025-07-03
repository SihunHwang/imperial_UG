#These functions read an ephemeris from the Horizons app and 
#decipher position data of the body at the given time step
#and output this data in a Nx4 list. Created by Debdut 
#Sengupta on 23 July 2022.

function loadEphemerides()

    #filenames=["Sun1.txt","Mercury1.txt","Venus1.txt","Earth1.txt","Mars1.txt",
    #"Jupiter1.txt","Saturn1.txt","Uranus1.txt","Neptune1.txt","Pluto1.txt","Ceres1.txt","Eris1.txt"]
    filenames = ["Sun2.txt", "Mercury2.txt", "Venus2.txt", "Earth2.txt", "Mars2.txt", "Jupiter2.txt", "Saturn2.txt", "Uranus2.txt", "Neptune2.txt"]

    ephemerides = [generateEphemeris(filename) for filename in filenames]

    return ephemerides
    
end

function findState(eph, t)

    #Simple identification, needs to be updated with Chebyshev Interpolation or similar

    #Can be speeded up - don't need to search through the ephemeris each time

    for row = 1:size(eph, 1)
        if (t < eph[row, 1]) #could replace t with (t + 1 day /2)
            return [eph[row, 2], eph[row, 3], eph[row, 4],
                eph[row, 5], eph[row, 6], eph[row, 7]]
        end
    end
    return [0.0, 0.0, 0.0]

end

function generateEphemeris(filename)
    #Ephemeris specifications: Runs from (latest) mission_start to (earliest)
    # mission_start+100. Should output T,X,Y,Z in the form of a Vector Table 
    # for the body in barycentric coordinates [SSB 500@0] with an interval of 
    # max 1 day and default settings.

    #filename type: "Earth#.txt" where # is replaced with trial number.

    file = open(filename, "r")

    global start = readline(file)
    while start != "\$\$SOE"
        start
        start = readline(file)
    end

    global row = 1
    global s = ' '
    global i = 0
    eph = zeros(100000, 7)
    s = readline(file)
    while s != "\$\$EOE"
        global T, X, Y, Z, VX, VY, VZ
        #Reading time
        i = 1
        T = ""
        while s[i] != ' '
            T = string(T, s[i])
            i = i + 1
        end
        s = readline(file)

        i = 5
        X = ""
        while s[i] != 'Y'
            X = string(X, s[i])
            i = i + 1
        end
        i = i + 3
        Y = ""
        while s[i] != 'Z'
            Y = string(Y, s[i])
            i = i + 1
        end
        i = i + 3
        Z = ""
        while i != length(s) + 1
            Z = string(Z, s[i])
            i = i + 1
        end

        s = readline(file)

        i = 5
        VX = ""
        while s[i] != 'V'
            VX = string(VX, s[i])
            i = i + 1
        end
        i = i + 3
        VY = ""
        while s[i] != 'V'
            VY = string(VY, s[i])
            i = i + 1
        end
        i = i + 3
        VZ = ""
        while i != length(s) + 1
            VZ = string(VZ, s[i])
            i = i + 1
        end

        #The ephemeride has the unit convention km-s although T is in h.
        eph[row, :] = [parse(Float64, T), parse(Float64, X) * 1000, parse(Float64, Y) * 1000,
            parse(Float64, Z) * 1000, parse(Float64, VX) * 1000, parse(Float64, VY) * 1000, parse(Float64, VZ) * 1000]
        row = row + 1

        readline(file)
        s = readline(file)
    end

    close(file)

    global j = 1
    while j <= size(eph, 1)
        if (eph[j, 1] == 0.0)
            break
        end
        global j = j + 1
    end

    eph = eph[1:j-1, :]
    #using Plots
    #plot!(eph[:,2],eph[:,3],eph[:,4])

    return eph

end
