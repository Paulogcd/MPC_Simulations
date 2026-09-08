using StatsFuns
using LinearAlgebra
using Distributions
using Random

function FnGridTrans(lx, ETrans1, ETrans2, RobustSDTrans1, RobustSDTrans2, pTrans1, ngpe, Twork)

    dTrans = MixtureModel(Normal[
            Normal(ETrans1, RobustSDTrans1),
            Normal(ETrans2, RobustSDTrans2)], [pTrans1, 1-pTrans1])

    Veps=var(dTrans)

    legrid = zeros(ngpe)
    ledist = zeros(ngpe)
    edist = zeros(Twork, ngpe)
    egrid = zeros(Twork, ngpe)

    legrid[1]=-lx*sqrt(Veps)
    legrid[ngpe] = lx*sqrt(Veps)
    lwidth = (legrid[ngpe]-legrid[1])/real(ngpe-1)
    ij=2
    while ij<=ngpe-1
        legrid[ij] = legrid[1] + lwidth*(ij-1)
        ij=ij+1
    end

    #! fill in probabilities using normal distribution
    ledist[1] = cdf(dTrans, (legrid[1]+0.5*lwidth)/sqrt(Veps))
    ij=2
    while ij<=ngpe-1
        ltemp1 = cdf(dTrans, (legrid[ij]+0.5*lwidth)/sqrt(Veps))
        ltemp3 = cdf(dTrans, (legrid[ij]-0.5*lwidth)/sqrt(Veps))
        ledist[ij] = ltemp1 - ltemp3
        ij=ij+1
    end
    ledist[ngpe]=1-cdf(dTrans, (legrid[ngpe]-0.5*lwidth)/sqrt(Veps)); # utiliser "upper()" pour plus de précision
    ledist = ledist ./ sum(ledist)

    #!find variance
    lvar = dot(legrid .^ 2, ledist) - dot(legrid, ledist) .^ 2

    #!moment
    FnGridTrans = (Veps - lvar)^2

    #!put values in globals
    it=1
    while it<=Twork
        edist[it, :] = ledist
        egrid[it, :] = legrid
        it=it+1
    end

    [FnGridTrans, edist, egrid]
end

export FnGridTrans