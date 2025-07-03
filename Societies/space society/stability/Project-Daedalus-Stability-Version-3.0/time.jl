#This file makes necessary time conversions for consistency of time 
# frames. It was created by Debdut Sengupta on 17th November 2022.


function yrsToSeconds(yrs)
    return yrs*31556926.0
end

function tToEphTime(t)
    mission_start = 2461041.5 #currently set to 2026-01-01. Needs to be updated
    return t / (3600 * 24) + mission_start
end