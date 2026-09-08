using Distributions
using Random

"""
Missing documentation
"""
function FnTaxParamNet(lstax, Twork, ngpz, ngpe, kappa, edist, egrid, zdist, zgrid, btax, ptax, pentax, popsize, theta, targetTaxToLabinc, Display)

    # RANDOM SEED

    rng = MersenneTwister(1234);

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

    g0=2.746
    g1=0.0624
    g2=-0.00167


    a=-2.495
    b=-0.1037
    c=-5.051
    d=-0.1087

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
            iep=1
            while iep<=ngpe
                g=g0+g1*it+g2*it*it
                xi=a+b*it+c*zgrid[it, iz]+d*it*zgrid[it, iz]
                p=exp(xi)/(1+exp(xi))
                if (p>=1)
                    dist = Binomial(1, 0.9)
                else
                    dist = Binomial(1, p)
                    display(p)
                end
                nu = rand(rng, dist, 1)[1]
                ygrid[it, iz, ie, iep] = (1-nu)*exp(g + egrid[it, ie] + theta*egrid[it, iep] + zgrid[it, iz])
                #=
                # get implied gross income at this point, to use in constructtion of soc sec system()
                lygross = ygrid[it,iz,ie,iep]/(1.0-pentax-btax)
                lygrossH = lygross*3.0
                lygrossL = 0.0
                lygrossacc = 1.0
                lygross=rtnewtGrossInc(lygross,ygrid[it,iz,ie,iep],lygrossL,lygrossH,lygrossacc,btax,ptax,stax,pentax)
                #rtnewtGrossInc[lnet,lxguess,lx1,lx2,xacc]
                #CALL rtnewtGrossInc[ygrid[it,iz,ie],lygross,lygrossL,lygrossH,lygrossacc]
                =#
                ypregrid[it, iz, ie, iep] = ygrid[it, iz, ie, iep]
                ltotlabincpre = ltotlabincpre + ypregrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it, iep]*popsize[it]
                ltotlabincpost = ltotlabincpost + ygrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it, iep]*popsize[it]
                FnTax = 0 #btax*(ypregrid[it,iz,ie,iep] - (ypregrid[it,iz,ie,iep]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it,iz,ie,iep]
                ltottax = 0 #ltottax + (FnTax - pentax*ypregrid[it,iz,ie,iep])*zdist[it,iz]*edist[it,ie]*edist[it,iep]*popsize[it]
                avearnspre[it] = avearnspre[it] + ypregrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avearnspre2[it] = avearnspre2[it] + (ypregrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avearnspost[it] = avearnspost[it] + ygrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avearnspost2[it] = avearnspost2[it] + (ygrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avlearnspre[it] = avlearnspre[it] + log(ypregrid[it, iz, ie, iep])*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avlearnspost[it] = avlearnspost[it] + log(ygrid[it, iz, ie, iep])*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                avlearnspost2[it] = avlearnspost2[it] + log(ygrid[it, iz, ie, iep]^2)*zdist[it, iz]*edist[it, ie]*edist[it, iep]
                iep=iep+1
            end
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
                    g=g0+g1*it+g2*it*it
                    xi=a+b*it+c*zgrid[it, iz]+d*it*zgrid[it, iz]
                    p=exp(xi)/(1+exp(xi))
                    if (p>=1)
                        dist = Binomial(1, 0.9)
                    else
                        dist = Binomial(1, p)
                    end
                    nu = rand(rng, dist, 1)[1]
                    ygrid[it, iz, ie, iep] = (1-nu)*exp(g + egrid[it, ie] + theta*egrid[it, iep] + zgrid[it, iz])
                    #=
                    # get implied gross income at this point, to use in constructtion of soc sec system()
                    lygross = ygrid[it,iz,ie,iep]/(1.0-pentax-btax)
                    lygrossH = lygross*3.0
                    lygrossL = 0.0
                    lygrossacc = 1.0
                    lygross=rtnewtGrossInc(lygross,ygrid[it,iz,ie,iep],lygrossL,lygrossH,lygrossacc,btax,ptax,stax,pentax)
                    #rtnewtGrossInc[lnet,lxguess,lx1,lx2,xacc]
                    #CALL rtnewtGrossInc[ygrid[it,iz,ie],lygross,lygrossL,lygrossH,lygrossacc]
                     =#
                    ypregrid[it, iz, ie, iep] = ygrid[it, iz, ie, iep]
                    ltotlabincpre = ltotlabincpre + ypregrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]*popsize[it]
                    ltotlabincpost = ltotlabincpost + ygrid[it, iz, ie, iep]*zdist[it, iz]*edist[it, ie]*edist[it-1, iep]*popsize[it]
                    FnTax = 0 #btax*(ypregrid[it,iz,ie,iep] - (ypregrid[it,iz,ie,iep]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it,iz,ie,iep]
                    ltottax = 0 #ltottax + (FnTax - pentax*ypregrid[it,iz,ie,iep])*zdist[it,iz]*edist[it,ie]*edist[it-1,iep]*popsize[it]
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

    FnTaxParamNet = ltottax/ltotlabincpre - targetTaxToLabinc
    totlabincpre = ltotlabincpre
    totlabincpost = ltotlabincpost

    if Display==1
        display(" Tax revenue / Pre-tax labor income: ")
        display((ltottax/ltotlabincpre)*100)
    end

    [FnTaxParamNet, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost];
end

export FnTaxParamNet