"""
Missing documentation.
"""
function Simulations(
    nsim, Ttot, rng, ngpal, aldist, Twork, edist)

    nnsim = zeros(nsim, Ttot);
    agesim = zeros(nsim, Ttot);
    
    # initial earnings
    
    alsimI[:] = sample(rng, collect(1:ngpal), weights(aldist[:]), nsim);
    zsimI[:, 1] = sample(rng, collect(1:ngpz), weights(zdist[1, :]), nsim);
    zsimIB = zeros(Int8, nsim, Twork);
    esimI[:, 1] = sample(rng, collect(1:ngpe), weights(edist[1, :]), nsim);
    # alsim = zeros(nsim);
    zsim = zeros(nsim, Twork);
    zsimB = zeros(nsim, Twork);
    esim = zeros(nsim, Twork);
    # alsim = algrid[alsimI[:]];
    zsim[:, 1] = zgrid[1, zsimI[:, 1]];
    esim[:, 1] = egrid[1, esimI[:, 1]];
    
    kappasim = zeros(nsim, Twork);
    kappasim[:, :] = ones(nsim, 1) * kappa';
    
    msimI = zeros(nsim, Twork);
    msimIB = zeros(nsim, Twork);
    ysim = zeros(nsim, Ttot);
    ysimB = zeros(nsim, Ttot);
    asimI = zeros(nsim, Ttot);
    asim = zeros(nsim, Ttot);
    asimB = zeros(nsim, Ttot);
    msim = zeros(nsim, Twork);
    msimB = zeros(nsim, Twork);
    
    trsim = zeros(nsim, Ttot);
    trsimB = zeros(nsim, Ttot);
    csim = zeros(nsim, Ttot);
    csimB = zeros(nsim, Ttot);
    csimT = zeros(nsim, Ttot);
    xsim = zeros(nsim, Ttot);
    xsimB = zeros(nsim, Ttot);
    
    ypresimB = ypresim;
    yavsimB = yavsim;
    
    Tbreak = zeros(nsim, Ttot);
    
    nu=0;
    
    # measurement error
    
    cmesim = randn!(rng, zeros(nsim, Ttot)) .* 2500;
    ymesim = randn!(rng, zeros(nsim, Ttot)) .* 2500;
    
    
    it_sim =1
    inn_sim=1
    while inn_sim<=nsim
        Tbreak[inn_sim, it_sim]=rand(rng, DiscreteUniform(2, Twork-6), 1)[1];
    
        # initial income
    
        nudist = Binomial(1, pzgrid[1, zsimI[inn_sim, 1]])
        nu = rand(rng, nudist, 1)[1]
        ypresim[inn_sim, 1] = ypregrid[1, zsimI[inn_sim, 1], esimI[inn_sim, 1], alsimI[inn_sim], nu+1]
        yavsim[inn_sim, 1] = min(ypresim[inn_sim, 1], pencap)
        lim, lp = FindLinProb1(mgrid[1, :], yavsim[inn_sim, 1])
        itemp = sample(rng, collect(1:2), weights(lp), 1)
        msimI[inn_sim, 1] = lim[itemp][1, 1]
        msim[inn_sim, 1] = mgrid[1, round(Int, msimI[inn_sim, 1])]
        ysim[inn_sim, 1] = ygrid[1, zsimI[inn_sim, 1], esimI[inn_sim, 1], alsimI[inn_sim], nu+1]
    
        # initial assets
    
        if RobustInitWealth==0
            itemp=sample(rng, collect(1:75), weights(initwealthdist[:, 2]), 1)
            asimI[inn_sim, 1] = itemp[1]
            asim[inn_sim, 1] = initwealthdist[itemp, 1][1]*ygrid[1, round(Int, (ngpz-1)/2-2), round(Int, (ngpe-1)/2-1), round(Int, (ngpal-1)/2-1)]
            asim[inn_sim, 1] = max(asim[inn_sim, 1], agrid[1, 1])
    
        elseif RobustInitWealth==1
            asim[:, 1] = 0.01*ones(nsim, 1); #0*ones(nsim,1); # -0.73*10030*ones(nsim,1) #-5*14030*ones(nsim,1)
        else
    
            itemp=sample(rng, collect(1:75), weights(initwealthdist[:, 2]), 1)
            asimI[inn_sim, 1] = itemp[1]
            asim[inn_sim, 1] = initwealthdist[itemp, 1][1]*ygrid[1, ngpz-2, 6, 2]
            asim[inn_sim, 1] = max(asim[inn_sim, 1], agrid[1, 1])
    
        end;
    
    
        nnsim[inn_sim, it_sim] = inn_sim;
        agesim[inn_sim, it_sim] = 1;
        if asim[inn_sim, it_sim]<0
            trsim[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
        else
            trsim[inn_sim, it_sim] = max(cfloor - (Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
        end
        csim[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimI[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asim[inn_sim, it_sim]+trsim[inn_sim, it_sim])
        csimT[inn_sim, it_sim] = (R-1) .* (asim[inn_sim, it_sim])+ysim[inn_sim, it_sim];
        if trsim[inn_sim, it_sim]>0.0
            csim[inn_sim, it_sim] = min(cfloor, csim[inn_sim, it_sim])
        end
        if asim[inn_sim, it_sim]<0
            xsim[inn_sim, it_sim] = Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] + trsim[inn_sim, it_sim]
        else
            xsim[inn_sim, it_sim] = Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] + trsim[inn_sim, it_sim]
        end
        asim[inn_sim, it_sim+1] = xsim[inn_sim, it_sim] - csim[inn_sim, it_sim]
        if asim[inn_sim, it_sim+1]<agrid[it_sim+1, 1]
            asim[inn_sim, it_sim+1] = agrid[it_sim+1, 1]
            csim[inn_sim, it_sim] = xsim[inn_sim, it_sim] - asim[inn_sim, it_sim+1]
        end
        if isnan(csim[inn_sim, it_sim])>0
            display("bad consumption")
        end
    
        inn_sim=inn_sim+1
    
    end;
    
    # working life
    
    psimI = zeros(nsim);
    psimIB=psimI;
    
    inn_sim=1
    while inn_sim<=nsim
    
        it_sim =2;
    
        while it_sim <=Tbreak[inn_sim, 1]-1
    
            Tbreak[inn_sim, it_sim]=Tbreak[inn_sim, 1];
            nnsim[inn_sim, it_sim] = inn_sim;
            agesim[inn_sim, it_sim] = it;
            zsimI[inn_sim, it_sim] = sample(rng, collect(1:ngpz), weights(ztrans[it-1, zsimI[inn_sim, it-1], :]), 1)[1, 1]
            esimI[inn_sim, it_sim] = sample(rng, collect(1:ngpe), weights(edist[it_sim, :]), 1)[1, 1]
            nudist = Binomial(1, pzgrid[it_sim, zsimI[inn_sim, it_sim]])
            nu = rand(rng, nudist, 1)[1]
            zsim[inn_sim, it_sim] = zgrid[it_sim, zsimI[inn_sim, it_sim]]
            esim[inn_sim, it_sim] = egrid[it_sim, esimI[inn_sim, it_sim]]
            ysim[inn_sim, it_sim] = ygrid[it_sim, zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1]
            ypresim[inn_sim, it_sim] = ypregrid[it_sim, zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1]
            yavsim[inn_sim, it_sim] = ((it-1)*msim[inn_sim, it-1] + min(ypresim[inn_sim, it_sim], pencap))/real(it)
            lim, lp = FindLinProb1(mgrid[it_sim, :], yavsim[inn_sim, it_sim])
            itemp = sample(rng, collect(1:2), weights(lp), 1)[1, 1]
            msimI[inn_sim, it_sim] = lim[itemp]
            msim[inn_sim, it_sim] = mgrid[it_sim, round(Int, msimI[inn_sim, it_sim])]
            if it_sim <Twork
                if asim[inn_sim, it_sim]<0
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
                else
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
                end
                csim[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimI[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asim[inn_sim, it_sim])
                csimT[inn_sim, it_sim] = (R-1) .* (asim[inn_sim, it_sim])+ysim[inn_sim, it_sim];
            elseif it_sim ==Twork
                psimI[inn_sim] = pind[round(Int, msimI[inn_sim, Twork]), zsimI[inn_sim, Twork], esimI[inn_sim, Twork], alsimI[inn_sim]]
                if asim[inn_sim, it_sim]<0
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agridret[1, 1, round(Int, psimI[inn_sim])]), 0.0)
                else
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agridret[1, 1, round(Int, psimI[inn_sim])]), 0.0)
                end
                csim[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimI[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asim[inn_sim, it_sim])
                csimT[inn_sim, it_sim] = (R-1) .* (asim[inn_sim, it_sim])+ysim[inn_sim, it_sim];
            end
            if asim[inn_sim, it_sim]<0
                xsim[inn_sim, it_sim] = Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] + trsim[inn_sim, it_sim]
            else
                xsim[inn_sim, it_sim] = Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] + trsim[inn_sim, it_sim]
            end
            asim[inn_sim, it_sim+1] = xsim[inn_sim, it_sim] - csim[inn_sim, it_sim]
    
            if it_sim <Twork
                if asim[inn_sim, it_sim+1]<agrid[it_sim+1, 1]
                    asim[inn_sim, it_sim+1] = agrid[it_sim+1, 1]
                    csim[inn_sim, it_sim] = xsim[inn_sim, it_sim] - asim[inn_sim, it_sim+1]
                end
            elseif it_sim ==Twork
                if asim[inn_sim, it_sim+1]<agridret[1, 1, round(Int, psimI[inn_sim])]
                    asim[inn_sim, it_sim+1] = agridret[1, 1, round(Int, psimI[inn_sim])]
                    csim[inn_sim, it_sim] = xsim[inn_sim, it_sim] - asim[inn_sim, it_sim+1]
                end
            end
    
            if isnan(csim[inn_sim, it_sim])>0
                display("bad consumption")
            end
    
            it_sim =it_sim+1
        end;
        asimB[:, it_sim]=asim[:, it_sim] .+ 0.5
    
        #edist and ztrans are what generate the trend
        while it_sim <=Twork
    
            Tbreak[inn_sim, it_sim]=Tbreak[inn_sim, 1];
            nnsim[inn_sim, it_sim] = inn_sim;
            agesim[inn_sim, it_sim] = it;
            zsimI[inn_sim, it_sim] = sample(rng, collect(1:ngpz), weights(ztrans[it-1, zsimI[inn_sim, it-1], :]), 1)[1, 1]
            esimI[inn_sim, it_sim] = sample(rng, collect(1:ngpe), weights(edist[it_sim, :]), 1)[1, 1]
            nudist = Binomial(1, pzgrid[it_sim, zsimI[inn_sim, it_sim]])
            nu = rand(rng, nudist, 1)[1]
            zsim[inn_sim, it_sim] = zgrid[it_sim, zsimI[inn_sim, it_sim]]
            esim[inn_sim, it_sim] = egrid[it_sim, esimI[inn_sim, it_sim]]
            ysim[inn_sim, it_sim] = ygrid[it_sim, zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1]
            ypresim[inn_sim, it_sim] = ypregrid[it_sim, zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1]
            yavsim[inn_sim, it_sim] = ((it-1)*msim[inn_sim, it-1] + min(ypresim[inn_sim, it_sim], pencap))/real(it)
            lim, lp = FindLinProb1(mgrid[it_sim, :], yavsim[inn_sim, it_sim])
            itemp = sample(rng, collect(1:2), weights(lp), 1)[1, 1]
            msimI[inn_sim, it_sim] = lim[itemp]
            msim[inn_sim, it_sim] = mgrid[it_sim, round(Int, msimI[inn_sim, it_sim])]
            if it_sim <Twork
                if asim[inn_sim, it_sim]<0
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
                else
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
                end
                csim[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimI[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asim[inn_sim, it_sim])
                csimT[inn_sim, it_sim] = (R-1) .* (asim[inn_sim, it_sim])+ysim[inn_sim, it_sim];
            elseif it_sim ==Twork
                psimI[inn_sim] = pind[round(Int, msimI[inn_sim, Twork]), zsimI[inn_sim, Twork], esimI[inn_sim, Twork], alsimI[inn_sim]]
                if asim[inn_sim, it_sim]<0
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agridret[1, 1, round(Int, psimI[inn_sim])]), 0.0)
                else
                    trsim[inn_sim, it_sim] = max(cfloor - (Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] - agridret[1, 1, round(Int, psimI[inn_sim])]), 0.0)
                end
                csim[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimI[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asim[inn_sim, it_sim])
                csimT[inn_sim, it_sim] = (R-1) .* (asim[inn_sim, it_sim])+ysim[inn_sim, it_sim];
            end
            if asim[inn_sim, it_sim]<0
                xsim[inn_sim, it_sim] = Rnetdebt*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] + trsim[inn_sim, it_sim]
            else
                xsim[inn_sim, it_sim] = Rnetsave*asim[inn_sim, it_sim] + ysim[inn_sim, it_sim] + trsim[inn_sim, it_sim]
            end
            asim[inn_sim, it_sim+1] = xsim[inn_sim, it_sim] - csim[inn_sim, it_sim]
    
            if it_sim <Twork
                if asim[inn_sim, it_sim+1]<agrid[it_sim+1, 1]
                    asim[inn_sim, it_sim+1] = agrid[it_sim+1, 1]
                    csim[inn_sim, it_sim] = xsim[inn_sim, it_sim] - asim[inn_sim, it_sim+1]
                end
            elseif it_sim ==Twork
                if asim[inn_sim, it_sim+1]<agridret[1, 1, round(Int, psimI[inn_sim])]
                    asim[inn_sim, it_sim+1] = agridret[1, 1, round(Int, psimI[inn_sim])]
                    csim[inn_sim, it_sim] = xsim[inn_sim, it_sim] - asim[inn_sim, it_sim+1]
                end
            end
    
            zsimB[inn_sim, it_sim] = zgrid[it_sim, zsimI[inn_sim, it_sim]]
            ysimB[inn_sim, it_sim] = ygrid[it_sim, zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1]
            ypresimB[inn_sim, it_sim] = ypregrid[it_sim, zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1]
            yavsimB[inn_sim, it_sim] = ((it-1)*msim[inn_sim, it-1] + min(ypresimB[inn_sim, it_sim], pencap))/real(it)
            lim, lp = FindLinProb1(mgrid[it_sim, :], yavsimB[inn_sim, it_sim])
            itemp = sample(rng, collect(1:2), weights(lp), 1)[1, 1]
            msimIB[inn_sim, it_sim] = lim[itemp]
            msimB[inn_sim, it_sim] = mgrid[it_sim, round(Int, msimI[inn_sim, it_sim])]
            if it_sim <Twork
                if asimB[inn_sim, it_sim]<0
                    trsimB[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asimB[inn_sim, it_sim] + ysimB[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
                else
                    trsimB[inn_sim, it_sim] = max(cfloor - (Rnetsave*asimB[inn_sim, it_sim] + ysimB[inn_sim, it_sim] - agrid[it_sim+1, 1]), 0.0)
                end
                csimB[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimIB[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asimB[inn_sim, it_sim])
            elseif it_sim ==Twork
                psimIB[inn_sim] = pind[round(Int, msimIB[inn_sim, Twork]), zsimI[inn_sim, Twork], esimI[inn_sim, Twork], alsimI[inn_sim]]
                if asimB[inn_sim, it_sim]<0
                    trsimB[inn_sim, it_sim] = max(cfloor - (Rnetdebt*asimB[inn_sim, it_sim] + ysimB[inn_sim, it_sim] - agridret[1, 1, round(Int, psimI[inn_sim])]), 0.0)
                else
                    trsimB[inn_sim, it_sim] = max(cfloor - (Rnetsave*asimB[inn_sim, it_sim] + ysimB[inn_sim, it_sim] - agridret[1, 1, round(Int, psimI[inn_sim])]), 0.0)
                end
                csimB[inn_sim, it_sim] = LinInterp1(agrid[it_sim, :], con[it_sim, :, round(Int, msimIB[inn_sim, it_sim]), zsimI[inn_sim, it_sim], esimI[inn_sim, it_sim], alsimI[inn_sim], nu+1], asimB[inn_sim, it_sim])
            end
            if asimB[inn_sim, it_sim]<0
                xsimB[inn_sim, it_sim] = Rnetdebt*asimB[inn_sim, it_sim] + ysimB[inn_sim, it_sim] + trsimB[inn_sim, it_sim]
            else
                xsimB[inn_sim, it_sim] = Rnetsave*asimB[inn_sim, it_sim] + ysimB[inn_sim, it_sim] + trsimB[inn_sim, it_sim]
            end
            asimB[inn_sim, it_sim+1] = xsimB[inn_sim, it_sim] - csimB[inn_sim, it_sim]
    
            if it_sim <Twork
                if asimB[inn_sim, it_sim+1]<agrid[it_sim+1, 1]
                    asimB[inn_sim, it_sim+1] = agrid[it_sim+1, 1]
                    csimB[inn_sim, it_sim] = xsimB[inn_sim, it_sim] - asimB[inn_sim, it_sim+1]
                end
            elseif it_sim ==Twork
                if asimB[inn_sim, it_sim+1]<agridret[1, 1, round(Int, psimIB[inn_sim])]
                    asimB[inn_sim, it_sim+1] = agridret[1, 1, round(Int, psimIB[inn_sim])]
                    csimB[inn_sim, it_sim] = xsimB[inn_sim, it_sim] - asimB[inn_sim, it_sim+1]
                end
            end
            # end if receiving shock
            if isnan(csim[inn_sim, it_sim])>0
                display("bad consumption")
            end
    
            it_sim =it_sim+1
        end
        inn_sim=inn_sim+1
    end
    
    # retirement
    
    it_sim = 1;
    
    while it_sim <= Tret
        inn_sim = 1
        while inn_sim <= nsim
            nnsim[inn_sim, Twork+it_sim] = inn_sim;
            agesim[inn_sim, Twork+it_sim] = Twork+it;
            ysim[inn_sim, Twork+it_sim] = pgrid[it_sim, round(Int, psimI[inn_sim])];
            ypresim[inn_sim, Twork+it_sim] = ppregrid[it_sim, round(Int, psimI[inn_sim])]
            if it_sim <Tret
                if asim[inn_sim, Twork+it_sim]<0
                    trsim[inn_sim, Twork+it_sim] = max(cfloor - (Rnetdebt*asim[inn_sim, Twork+it_sim] + ysim[inn_sim, Twork+it_sim] - agridret[it_sim+1, 1, round(Int, psimI[inn_sim])]), 0.0)
                else
                    trsim[inn_sim, Twork+it_sim] = max(cfloor - (Rnetsave*asim[inn_sim, Twork+it_sim] + ysim[inn_sim, Twork+it_sim] - agridret[it_sim+1, 1, round(Int, psimI[inn_sim])]), 0.0)
                end
                csim[inn_sim, Twork+it_sim] = LinInterp1(agridret[it_sim, :, round(Int, psimI[inn_sim])], conret[it_sim, :, round(Int, psimI[inn_sim])], asim[inn_sim, Twork+it_sim])
                csimT[inn_sim, Twork+it_sim] = (R-1) .* (asim[inn_sim, Twork+it_sim])+ysim[inn_sim, Twork+it_sim];
                if asim[inn_sim, Twork+it_sim]<0
                    xsim[inn_sim, Twork+it_sim] = Rnetdebt*asim[inn_sim, Twork+it_sim] + ysim[inn_sim, Twork+it_sim] + trsim[inn_sim, Twork+it_sim]
                else
                    xsim[inn_sim, Twork+it_sim] = Rnetsave*asim[inn_sim, Twork+it_sim] + ysim[inn_sim, Twork+it_sim] + trsim[inn_sim, Twork+it_sim]
                end
                asim[inn_sim, Twork+it_sim+1] = (xsim[inn_sim, Twork+it_sim] - csim[inn_sim, Twork+it_sim])/annprem[it_sim]
            end
            if it_sim <Tret
                if asim[inn_sim, Twork+it_sim+1]<agridret[it_sim+1, 1, round(Int, psimI[inn_sim])]
                    asim[inn_sim, Twork+it_sim+1] = agridret[it_sim+1, 1, round(Int, psimI[inn_sim])]
                    csim[inn_sim, Twork+it_sim] = xsim[inn_sim, Twork+it_sim] - asim[inn_sim, Twork+it_sim+1]*annprem[it_sim]
                end
            end
            inn_sim=inn_sim+1
        end
        it_sim =it_sim+1
    end;
    
    csim_withme=csim .* exp.(randn!(rng, zeros(nsim, Ttot)) .* 0.5) #.+cmesim;
    ysim_withme=ysim .* exp.(randn!(rng, zeros(nsim, Ttot)) .* 0.5) #.+ymesim;
    
    alsim = zeros(nsim, Ttot);
    alsim = algrid[alsimI[:]]*ones(Ttot)';
    
    # Mean values for all simulations
    csim_mean = zeros(Ttot);
    csimB_mean = zeros(Ttot);
    csimT_mean = zeros(Ttot);
    asim_mean = zeros(Ttot);
    ysim_mean = zeros(Ttot);
    ysimB_mean = zeros(Ttot);
    
    for age=1:(Ttot-1)
        csim_mean[age] = mean(csim[:, age]);
        csimB_mean[age] = mean(csimB[:, age]);
        csimT_mean[age] = mean(csimT[:, age]);
        asim_mean[age] = mean(asim[:, age]);
        ysim_mean[age] = mean(ysim[:, age]);
        ysimB_mean[age] = mean(ysimB[:, age]);
        print(" \n")
        print("Constraints? min - max\n")
        print(minimum(asim[:, age]));
        print(" \n")
        print(maximum(asim[:, age]));
        print(" \n")
    end
end;
export Simulations;