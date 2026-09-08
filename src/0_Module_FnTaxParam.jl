"""
Missing documentation
"""
function FnTaxParam(lstax, Twork, ngpz, ngpe, kappa, edist, egrid, zdist, zgrid, btax, ptax, pentax, popsize, theta, targetTaxToLabinc, Display)

    avearnspre = zeros(Twork)
    avearnspre2 = zeros(Twork)
    avearnspost = zeros(Twork)
    avearnspost2 = zeros(Twork)
    avlearnspre = zeros(Twork)
    avlearnspre2 = zeros(Twork)
    avlearnspost = zeros(Twork)
    avlearnspost2 = zeros(Twork)
    ygrid = zeros(Twork, ngpz, ngpe, ngpe)
    ypregrid = zeros(Twork, ngpz, ngpe, ngpe)

    stax = lstax

    ltottax = 0.0
    ltotlabincpre = 0.0
    ltotlabincpost = 0.0

    it=1
    avearnspre[it] = 0.0
    avearnspre2[it] = 0.0
    avearnspost[it] = 0.0
    avearnspost2[it] = 0.0
    avlearnspre[it] = 0.0
    avlearnspre2[it] = 0.0
    avlearnspost[it] = 0.0
    avlearnspost2[it] = 0.0
    iz=1
    while iz<=ngpz
        ie=1
        while ie<=ngpe
            ypregrid[it, iz, ie, ie] = exp(kappa[it] + egrid[it, ie] + theta*egrid[it, ie] + zgrid[it, iz])
            FnTax = btax*(ypregrid[it, iz, ie, ie] - (ypregrid[it, iz, ie, ie]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it, iz, ie, ie]
            ygrid[it, iz, ie, ie] = ypregrid[it, iz, ie, ie] - FnTax

            ltotlabincpre = ltotlabincpre + ypregrid[it, iz, ie, ie]*zdist[it, iz]*edist[it, ie]*popsize[it]
            ltotlabincpost = ltotlabincpost + ygrid[it, iz, ie, ie]*zdist[it, iz]*edist[it, ie]*popsize[it]

            ltottax = ltottax + (FnTax - pentax*ypregrid[it, iz, ie, ie])*zdist[it, iz]*edist[it, ie]*popsize[it]
            avearnspre[it] = avearnspre[it] + ypregrid[it, iz, ie, ie]*zdist[it, iz]*edist[it, ie]
            avearnspre2[it] = avearnspre2[it] + (ypregrid[it, iz, ie, ie]^2)*zdist[it, iz]*edist[it, ie]
            avearnspost[it] = avearnspost[it] + ygrid[it, iz, ie, ie]*zdist[it, iz]*edist[it, ie]
            avearnspost2[it] = avearnspost2[it] + (ygrid[it, iz, ie, ie]^2)*zdist[it, iz]*edist[it, ie]
            avlearnspre[it] = avlearnspre[it] + log(ypregrid[it, iz, ie, ie])*zdist[it, iz]*edist[it, ie]
            avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it, iz, ie, ie]^2)*zdist[it, iz]*edist[it, ie]
            avlearnspost[it] = avlearnspost[it] + log(ygrid[it, iz, ie, ie])*zdist[it, iz]*edist[it, ie]
            avlearnspost2[it] = avlearnspost2[it] + log(ygrid[it, iz, ie, ie]^2)*zdist[it, iz]*edist[it, ie]
            ie=ie+1
        end
        iz=iz+1
    end

    it=2
    while it<=Twork
        avearnspre[it] = 0.0
        avearnspre2[it] = 0.0
        avearnspost[it] = 0.0
        avearnspost2[it] = 0.0
        avlearnspre[it] = 0.0
        avlearnspre2[it] = 0.0
        avlearnspost[it] = 0.0
        avlearnspost2[it] = 0.0
        iz=1
        while iz<=ngpz
            ie=1
            while ie<=ngpe
                iep=1
                while iep<=ngpe
                    ypregrid[it, iz, ie, iep] = exp(kappa[it] + egrid[it, ie] + theta*egrid[it-1, iep] + zgrid[it, iz])
                    FnTax = btax*(ypregrid[it, iz, ie, iep] - (ypregrid[it, iz, ie, iep]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it, iz, ie, iep]
                    ygrid[it, iz, ie, iep] = ypregrid[it, iz, ie, iep] - FnTax

                    ltotlabincpre = ltotlabincpre + ypregrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]*popsize[it]
                    ltotlabincpost = ltotlabincpost + ygrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]*popsize[it]

                    ltottax = ltottax + (FnTax - pentax*ypregrid[it, iz, ie, iep])*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]*popsize[it]
                    avearnspre[it] = avearnspre[it] + ypregrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avearnspre2[it] = avearnspre2[it] + (ypregrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avearnspost[it] = avearnspost[it] + ygrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avearnspost2[it] = avearnspost2[it] + (ygrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avlearnspre[it] = avlearnspre[it] + log(ypregrid[it, iz, ie, iep])*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avlearnspost[it] = avlearnspost[it] + log(ygrid[it, iz, ie, iep])*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    avlearnspost2[it] = avlearnspost2[it] + log(ygrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]
                    iep=iep+1
                end
                ie=ie+1
            end
            iz=iz+1
        end
        it=it+1
    end

    varearnspre = avearnspre2 - avearnspre .^ 2
    varearnspost = avearnspost2 - avearnspost .^ 2
    varlearnspre = avlearnspre2 - avlearnspre .^ 2
    varlearnspost = avlearnspost2 - avlearnspost .^ 2

    FnTaxParam = ltottax/ltotlabincpre - targetTaxToLabinc
    totlabincpre = ltotlabincpre
    totlabincpost = ltotlabincpost

    if Display==1
        display(" Tax revenue / Pre-tax labor income: ")
        display((ltottax/ltotlabincpre)*100)
    end

    [FnTaxParam, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost];
end

export FnTaxParam