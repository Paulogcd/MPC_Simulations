using StatsFuns
using LinearAlgebra
using Distributions
using Random

"""
Missing documentation
"""
function FnGridKEntrep(lx, SDAlpha, ngpal)

    dAlpha = Normal(0, SDAlpha)

    legrid = zeros(ngpal)
    ledist = zeros(ngpal)
    aldist = zeros(ngpal)
    algrid = zeros(ngpal)

    legrid[1]=-lx * SDAlpha
    legrid[ngpal] = lx * SDAlpha
    lwidth = (legrid[ngpal]-legrid[1])/real(ngpal-1)
    ij=2
    while ij<=ngpal-1
        legrid[ij] = legrid[1] + lwidth*(ij-1)
        ij=ij+1
    end

    #! fill in probabilities using normal distribution
    ledist[1] = cdf(dAlpha, (legrid[1]+0.5*lwidth) / SDAlpha)
    ij=2
    while ij<=ngpal-1
        ltemp1 = cdf(dAlpha, (legrid[ij]+0.5*lwidth) / SDAlpha)
        ltemp3 = cdf(dAlpha, (legrid[ij]-0.5*lwidth) / SDAlpha)
        ledist[ij] = ltemp1 - ltemp3
        ij=ij+1
    end
    ledist[ngpal]=1-cdf(dAlpha, (legrid[ngpal]-0.5*lwidth)/SDAlpha); # utiliser "upper()" pour plus de précision
    ledist = ledist ./ sum(ledist)

    #!find variance
    lvar = dot(legrid .^ 2, ledist) - dot(legrid, ledist) .^ 2

    #!moment
    FnGridAlpha = (SDAlpha^2 - lvar)^2

    #!put values in globals
    aldist[:] = ledist
    algrid[:] = legrid


    [FnGridAlpha, aldist, algrid]
end
