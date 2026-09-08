using StatsFuns
using LinearAlgebra
using Distributions
using Random

""" 
Missing documentation.
"""
function FnGridPerm(lx, EPerm1, EPerm2, RobustSDPerm1, RobustSDPerm2, pPerm1, Twork, varz, ngpz, rho, SDz0)

    dPerm = MixtureModel(Normal[
            Normal(EPerm1, RobustSDPerm1),
            Normal(EPerm2, RobustSDPerm2)], [pPerm1, 1-pPerm1])
    Vetavec=var(dPerm)*ones(Twork-1)

    zgrid = zeros(Twork, ngpz)
    lwidth = zeros(Twork)
    ztrans = zeros(Twork-1, ngpz, ngpz)
    zdist = zeros(Twork, ngpz)
    lvar = zeros(Twork)

    #! get boundaries and fill in with equally spaced points
    it=1
    while it<=Twork
        zgrid[it, 1] = -lx*sqrt(varz[it])
        zgrid[it, ngpz] = lx*sqrt(varz[it])
        lwidth[it] = (zgrid[it, ngpz]-zgrid[it, 1])/real(ngpz-1)
        iz1=2
        while iz1<=ngpz-1
            zgrid[it, iz1] = zgrid[it, 1] + lwidth[it]*(iz1-1)
            iz1=iz1+1
        end
        it=it+1
    end

    #! fill in transition matrix using normal distribution
    it=1
    while it<=Twork-1
        iz1=1
        while iz1<=ngpz
            ztrans[it, iz1, 1]=cdf(dPerm, (zgrid[it+1, 1]+0.5*lwidth[it+1]-rho*zgrid[it, iz1])/sqrt(Vetavec[it]))
            iz2=2
            while iz2<=ngpz
                ltemp1=cdf(dPerm, (zgrid[it+1, iz2]+0.5*lwidth[it+1]-rho*zgrid[it, iz1])/sqrt(Vetavec[it]))
                ltemp3=cdf(dPerm, (zgrid[it+1, iz2]-0.5*lwidth[it+1]-rho*zgrid[it, iz1])/sqrt(Vetavec[it]))
                ztrans[it, iz1, iz2] = ltemp1 - ltemp3
                iz2=iz2+1
            end
            ztrans[it, iz1, ngpz]=1-cdf(dPerm, (zgrid[it+1, ngpz]-0.5*lwidth[it+1]-rho*zgrid[it, iz1])/sqrt(Vetavec[it]))
            ztrans[it, iz1, :] = ztrans[it, iz1, :] ./ sum(ztrans[it, iz1, :])
            iz1=iz1+1
        end
        it=it+1
    end

    #!find distribution at first period
    zdist[1, 1]=normcdf((zgrid[1, 1]+0.5*lwidth[1])/SDz0)
    iz1=2
    while iz1<=ngpz-1
        ltemp1=normcdf((zgrid[1, iz1]+0.5*lwidth[1])/SDz0)
        ltemp3=normcdf((zgrid[1, iz1]-0.5*lwidth[1])/SDz0)
        zdist[1, iz1] = ltemp1 - ltemp3
        iz1=iz1+1
    end
    zdist[1, ngpz]=1-normcdf((zgrid[1, ngpz]-0.5*lwidth[1])/SDz0)
    zdist[1, :] = zdist[1, :] ./ sum(zdist[1, :])

    #!find unconditional distributions
    it=2

    while it<=Twork
        A = dropdims(ztrans[it-1, :, :]; dims=tuple(findall(size(ztrans[it-1, :, :]) .== 1)...))
        zdist[it, :] = transpose(zdist[it-1, :])*A
        zdist[it, :] = zdist[it, :] ./ sum(zdist[it, :])
        it=it+1
    end

    #!find variance
    it=1
    while it<=Twork
        lvar[it] = dot(zgrid[it, :] .^ 2, zdist[it, :]) - dot(zgrid[it, :], zdist[it, :]) .^ 2
        it=it+1
    end

    #!moment
    FnGridPerm = sum((varz - lvar) .^ 2)

    #!put values in globals
    global varzapprox = lvar

    [FnGridPerm, varzapprox, zdist, zgrid, ztrans]
end

export FnGridPerm