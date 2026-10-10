using StatsBase
export Simulations

"""
Missing documentation.
"""
function Simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)
    
    Initialisation_simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid);

    # Initial state (do t=1 separately)
    Initial_state_simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid);

    # working life (15)
    Working_life_simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid);
                                                        
    # Retirement
    Retirement_simulation(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid);

    # Return formatted result
    results = Result_simulation(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid);

    @info ("\nSimulations done.")

    return results
end

function Initialisation_simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)

    global nnsim  = zeros(nsim, Ttot);
    global agesim = zeros(nsim, Ttot);
    
    # initial earnings

    global alsimI[:]=sample(rng, collect(1:ngpal), weights(aldist[:]),nsim);
    global zsimI[:,1]=sample(rng, collect(1:ngpz), weights(zdist[1,:]), nsim);
    global zsimIB= zeros(Int8, nsim, Twork);
    global esimI[:,1]=sample(rng, collect(1:ngpe), weights(edist[1,:]), nsim);
    # global #alsim = zeros(nsim);
    global zsim = zeros(nsim, Twork);
    global zsimB = zeros(nsim, Twork);
    global esim = zeros(nsim, Twork);
    # global #alsim = algrid[alsimI[:]];
    global zsim[:,1] = zgrid[1,zsimI[:,1]];
    global esim[:,1] = egrid[1,esimI[:,1]];

    global kappasim=zeros(nsim,Twork);
    global kappasim[:,:] = ones(nsim,1)*kappa';

    global msimI   = zeros(nsim, Twork);
    global msimIB  = zeros(nsim, Twork);
    global ysim    = zeros(nsim, Ttot);
    global ysimB   = zeros(nsim, Ttot);
    global asimI   = zeros(nsim, Ttot);
    global asim    = zeros(nsim, Ttot);
    global asimB   = zeros(nsim, Ttot);
    global msim    = zeros(nsim, Twork);
    global msimB   = zeros(nsim, Twork);

    global trsim       = zeros(nsim,Ttot);
    global trsimB      = zeros(nsim,Ttot);
    global csim        = zeros(nsim, Ttot);
    global csimB       = zeros(nsim, Ttot);
    global csimT       = zeros(nsim,Ttot);
    global xsim        = zeros(nsim, Ttot);
    global xsimB       = zeros(nsim, Ttot);
    global ypresimB    = ypresim;
    global yavsimB     = yavsim;
    global Tbreak      = zeros(nsim, Ttot);
    global nu          = 0;
end

"""

    Initial_state_simulations()


    Helper function 2/5 for the simulations.
    
    This function computes the simulated variables at the initial state of the agents.
"""
function Initial_state_simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)
    
    it=1

    Threads.@threads for inn in 1:nsim
        
        rng = MersenneTwister(1234 + inn)
        Tbreak[inn,it] = rand(rng, DiscreteUniform(2,Twork-6))
        
        # initial income

        nudist          = Binomial(1,pzgrid[1,zsimI[inn,1]])
        nu              = rand(rng, nudist)
        ypresim[inn,1]  = ypregrid[1,zsimI[inn,1],esimI[inn,1],alsimI[inn],nu+1]
        yavsim[inn,1]   = min(ypresim[inn,1],pencap)
        lim,lp          = FindLinProb1(mgrid[1,:],yavsim[inn,1])
        itemp           = sample(rng, collect(1:2), weights(lp))
        msimI[inn,1]    = lim[itemp][1,1]
        msim[inn,1]     = mgrid[1,round(Int,msimI[inn,1])]
        ysim[inn,1]     = ygrid[1,zsimI[inn,1],esimI[inn,1],alsimI[inn],nu+1]

        # initial assets

        if RobustInitWealth == 0
            itemp=sample(rng, collect(1:75), weights(initwealthdist[:,2]))
            asimI[inn,1]    = itemp[1]
            asim[inn,1]     = initwealthdist[itemp,1][1]*ygrid[1,round(Int,(ngpz-1)/2-2),round(Int,(ngpe-1)/2-1),round(Int,(ngpal-1)/2-1)]
            asim[inn,1]     = max(asim[inn,1], agrid[1,1])

        elseif RobustInitWealth == 1
            asim[:,1] = 0.01 *ones(nsim,1); #0*ones(nsim,1); # -0.73*10030*ones(nsim,1) #-5*14030*ones(nsim,1)
        else
            itemp=sample(rng, collect(1:75), weights(initwealthdist[:,2]))
            asimI[inn,1] = itemp[1]
            asim[inn,1] = initwealthdist[itemp,1][1]*ygrid[1,ngpz-2,6,2]
            asim[inn,1] = max(asim[inn,1], agrid[1,1])
        end;

        nnsim[inn,it] = inn;
        agesim[inn,it] = 1;
        if asim[inn,it]<0
            trsim[inn,it]    = max( cfloor - (Rnetdebt*asim[inn,it]+ ysim[inn,it] - agrid[it+1,1]) ,0.0)
        else
            trsim[inn,it]    = max( cfloor - (Rnetsave*asim[inn,it]+ ysim[inn,it] - agrid[it+1,1]) ,0.0)
        end
        csim[inn,it] = LinInterp1(agrid[it,:],con[it,:,round(Int,msimI[inn,it]),zsimI[inn,it],esimI[inn,it],alsimI[inn],nu+1],asim[inn,it]+trsim[inn,it])
        csimT[inn,it] = (R-1).*(asim[inn,it])+ysim[inn,it];
        if trsim[inn,it]>0.0
            csim[inn,it] = min(cfloor, csim[inn,it])
        end
        if asim[inn,it]<0
            xsim[inn,it]     = Rnetdebt*asim[inn,it]+ ysim[inn,it] + trsim[inn,it]
        else
            xsim[inn,it]     = Rnetsave*asim[inn,it]+ ysim[inn,it] + trsim[inn,it]
        end
        asim[inn,it+1]   = xsim[inn,it] - csim[inn,it]
        if asim[inn,it+1]<agrid[it+1,1]
            asim[inn,it+1] = agrid[it+1,1]
            csim[inn,it] = xsim[inn,it] - asim[inn,it+1]
        end
        if isnan(csim[inn,it])>0
            display("bad consumption")
        end
    end
end

function Working_life_simulations(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)
    
    global psimI = zeros(nsim);
    global psimIB=psimI;

    Threads.@threads for inn in 1:nsim

        rng = MersenneTwister(1234 + inn)
        
        it_3=2;

        while it_3<=Tbreak[inn,1]-1
            
            Tbreak[inn,it_3]=Tbreak[inn,1];
            nnsim[inn,it_3] = inn;
            agesim[inn,it_3]= it_3;
            zsimI[inn,it_3]= sample(rng, collect(1:ngpz), weights(ztrans[it_3-1,zsimI[inn,it_3-1],:]), 1)[1,1]
            esimI[inn,it_3]= sample(rng, collect(1:ngpe), weights(edist[it_3,:]), 1)[1,1]
            nudist = Binomial(1,pzgrid[it_3,zsimI[inn,it_3]])
            nu = rand(rng, nudist, 1)[1]
            zsim[inn,it_3] = zgrid[it_3,zsimI[inn,it_3]]
            esim[inn,it_3] = egrid[it_3,esimI[inn,it_3]]
            ysim[inn,it_3] = ygrid[it_3,zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1]
            ypresim[inn,it_3] = ypregrid[it_3,zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1]
            yavsim[inn,it_3] = ((it_3-1)*msim[inn,it_3-1] + min(ypresim[inn,it_3],pencap))/real(it_3)
            lim,lp = FindLinProb1(mgrid[it_3,:],yavsim[inn,it_3])
            it_3emp= sample(rng, collect(1:2), weights(lp), 1)[1,1]
            msimI[inn,it_3] = lim[it_3emp]
            msim[inn,it_3] = mgrid[it_3,round(Int,msimI[inn,it_3])]
            if it_3<Twork
                if asim[inn,it_3]<0
                    trsim[inn,it_3]    = max( cfloor - (Rnetdebt*asim[inn,it_3]+ ysim[inn,it_3] - agrid[it_3+1,1]) ,0.0)
                else
                    trsim[inn,it_3]    = max( cfloor - (Rnetsave*asim[inn,it_3]+ ysim[inn,it_3] - agrid[it_3+1,1]) ,0.0)
                end
                csim[inn,it_3] = LinInterp1(agrid[it_3,:],con[it_3,:,round(Int,msimI[inn,it_3]),zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1],asim[inn,it_3])
                csimT[inn,it_3] = (R-1).*(asim[inn,it_3])+ysim[inn,it_3];
            elseif it_3==Twork
                psimI[inn] = pind[round(Int, msimI[inn,Twork]),zsimI[inn,Twork],esimI[inn,Twork],alsimI[inn]]
                if asim[inn,it_3]<0
                    trsim[inn,it_3]    = max( cfloor - (Rnetdebt*asim[inn,it_3]+ ysim[inn,it_3] - agridret[1,1,round(Int,psimI[inn])] ) ,0.0)
                else
                    trsim[inn,it_3]    = max( cfloor - (Rnetsave*asim[inn,it_3]+ ysim[inn,it_3] - agridret[1,1,round(Int,psimI[inn])] ) ,0.0)
                end
                csim[inn,it_3] = LinInterp1(agrid[it_3,:],con[it_3,:,round(Int, msimI[inn,it_3]),zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1],asim[inn,it_3])
                csimT[inn,it_3] = (R-1).*(asim[inn,it_3])+ysim[inn,it_3];
            end
            if asim[inn,it_3]<0
                xsim[inn,it_3]     = Rnetdebt*asim[inn,it_3]+ ysim[inn,it_3] + trsim[inn,it_3]
            else
                xsim[inn,it_3]     = Rnetsave*asim[inn,it_3]+ ysim[inn,it_3] + trsim[inn,it_3]
            end
            asim[inn,it_3+1]   = xsim[inn,it_3] - csim[inn,it_3]
            
            if it_3<Twork
                if asim[inn,it_3+1]<agrid[it_3+1,1]
                    asim[inn,it_3+1] = agrid[it_3+1,1]
                    csim[inn,it_3] = xsim[inn,it_3] - asim[inn,it_3+1]
                end
            elseif it_3==Twork
                if asim[inn,it_3+1]<agridret[1,1,round(Int,psimI[inn])]
                    asim[inn,it_3+1] = agridret[1,1,round(Int,psimI[inn])]
                    csim[inn,it_3] = xsim[inn,it_3] - asim[inn,it_3+1]
                end
            end
            
            if isnan(csim[inn,it_3])>0
                display("bad consumption")
            end

            it_3=it_3+1
        
            asimB[:,it_3]=asim[:,it_3].+0.5 #=1*ysim[:,it_3-1]=#
        end

        #edist and ztrans are what generate the trend
        while it_3<=Twork
            
            Tbreak[inn,it_3]=Tbreak[inn,1];
            nnsim[inn,it_3] = inn;
            agesim[inn,it_3]= it_3;
            zsimI[inn,it_3]= sample(rng, collect(1:ngpz), weights(ztrans[it_3-1,zsimI[inn,it_3-1],:]), 1)[1,1]
            esimI[inn,it_3]= sample(rng, collect(1:ngpe), weights(edist[it_3,:]), 1)[1,1]
            nudist = Binomial(1,pzgrid[it_3,zsimI[inn,it_3]])
            nu = rand(rng, nudist, 1)[1]
            zsim[inn,it_3] = zgrid[it_3,zsimI[inn,it_3]]
            esim[inn,it_3] = egrid[it_3,esimI[inn,it_3]]
            ysim[inn,it_3] = ygrid[it_3,zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1]
            ypresim[inn,it_3] = ypregrid[it_3,zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1]
            yavsim[inn,it_3] = ((it_3-1)*msim[inn,it_3-1] + min(ypresim[inn,it_3],pencap))/real(it_3)
            lim,lp = FindLinProb1(mgrid[it_3,:],yavsim[inn,it_3])
            it_3emp= sample(rng, collect(1:2), weights(lp), 1)[1,1]
            msimI[inn,it_3] = lim[it_3emp]
            msim[inn,it_3] = mgrid[it_3,round(Int,msimI[inn,it_3])]
            if it_3<Twork
                if asim[inn,it_3]<0
                    trsim[inn,it_3]    = max( cfloor - (Rnetdebt*asim[inn,it_3]+ ysim[inn,it_3] - agrid[it_3+1,1]) ,0.0)
                else
                    trsim[inn,it_3]    = max( cfloor - (Rnetsave*asim[inn,it_3]+ ysim[inn,it_3] - agrid[it_3+1,1]) ,0.0)
                end
                csim[inn,it_3] = LinInterp1(agrid[it_3,:],con[it_3,:,round(Int,msimI[inn,it_3]),zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1],asim[inn,it_3])
                csimT[inn,it_3] = (R-1).*(asim[inn,it_3])+ysim[inn,it_3];
            elseif it_3==Twork
                psimI[inn] = pind[round(Int, msimI[inn,Twork]),zsimI[inn,Twork],esimI[inn,Twork],alsimI[inn]]
                if asim[inn,it_3]<0
                    trsim[inn,it_3]    = max( cfloor - (Rnetdebt*asim[inn,it_3]+ ysim[inn,it_3] - agridret[1,1,round(Int,psimI[inn])] ) ,0.0)
                else
                    trsim[inn,it_3]    = max( cfloor - (Rnetsave*asim[inn,it_3]+ ysim[inn,it_3] - agridret[1,1,round(Int,psimI[inn])] ) ,0.0)
                end
                csim[inn,it_3] = LinInterp1(agrid[it_3,:],con[it_3,:,round(Int, msimI[inn,it_3]),zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1],asim[inn,it_3])
                csimT[inn,it_3] = (R-1).*(asim[inn,it_3])+ysim[inn,it_3];
            end
            if asim[inn,it_3]<0
                xsim[inn,it_3]     = Rnetdebt*asim[inn,it_3]+ ysim[inn,it_3] + trsim[inn,it_3]
            else
                xsim[inn,it_3]     = Rnetsave*asim[inn,it_3]+ ysim[inn,it_3] + trsim[inn,it_3]
            end
            asim[inn,it_3+1]   = xsim[inn,it_3] - csim[inn,it_3]
            
            if it_3<Twork
                if asim[inn,it_3+1]<agrid[it_3+1,1]
                    asim[inn,it_3+1] = agrid[it_3+1,1]
                    csim[inn,it_3] = xsim[inn,it_3] - asim[inn,it_3+1]
                end
            elseif it_3==Twork
                if asim[inn,it_3+1]<agridret[1,1,round(Int,psimI[inn])]
                    asim[inn,it_3+1] = agridret[1,1,round(Int,psimI[inn])]
                    csim[inn,it_3] = xsim[inn,it_3] - asim[inn,it_3+1]
                end
            end
            
            zsimB[inn,it_3] = zgrid[it_3,zsimI[inn,it_3]]
            ysimB[inn,it_3] = ygrid[it_3,zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1]
            ypresimB[inn,it_3] = ypregrid[it_3,zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1]
            yavsimB[inn,it_3] = ((it_3-1)*msim[inn,it_3-1] + min(ypresimB[inn,it_3],pencap))/real(it_3)
            lim,lp = FindLinProb1(mgrid[it_3,:],yavsimB[inn,it_3])
            it_3emp= sample(rng, collect(1:2), weights(lp), 1)[1,1]
            msimIB[inn,it_3] = lim[it_3emp]
            msimB[inn,it_3] = mgrid[it_3,round(Int,msimI[inn,it_3])]
            if it_3<Twork
                if asimB[inn,it_3]<0
                    trsimB[inn,it_3]    = max( cfloor - (Rnetdebt*asimB[inn,it_3]+ ysimB[inn,it_3] - agrid[it_3+1,1]) ,0.0)
                else
                    trsimB[inn,it_3]    = max( cfloor - (Rnetsave*asimB[inn,it_3]+ ysimB[inn,it_3] - agrid[it_3+1,1]) ,0.0)
                end
                csimB[inn,it_3] = LinInterp1(agrid[it_3,:],con[it_3,:,round(Int,msimIB[inn,it_3]),zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1],asimB[inn,it_3])
            elseif it_3==Twork
                psimIB[inn] = pind[round(Int, msimIB[inn,Twork]),zsimI[inn,Twork],esimI[inn,Twork],alsimI[inn]]
                if asimB[inn,it_3]<0
                    trsimB[inn,it_3]    = max( cfloor - (Rnetdebt*asimB[inn,it_3]+ ysimB[inn,it_3] - agridret[1,1,round(Int,psimI[inn])] ) ,0.0)
                else
                    trsimB[inn,it_3]    = max( cfloor - (Rnetsave*asimB[inn,it_3]+ ysimB[inn,it_3] - agridret[1,1,round(Int,psimI[inn])] ) ,0.0)
                end
                csimB[inn,it_3] = LinInterp1(agrid[it_3,:],con[it_3,:,round(Int, msimIB[inn,it_3]),zsimI[inn,it_3],esimI[inn,it_3],alsimI[inn],nu+1],asimB[inn,it_3])
            end
            if asimB[inn,it_3]<0
                xsimB[inn,it_3]     = Rnetdebt*asimB[inn,it_3]+ ysimB[inn,it_3] + trsimB[inn,it_3]
            else
                xsimB[inn,it_3]     = Rnetsave*asimB[inn,it_3]+ ysimB[inn,it_3] + trsimB[inn,it_3]
            end
            asimB[inn,it_3+1]   = xsimB[inn,it_3] - csimB[inn,it_3]
            
            if it_3<Twork
                if asimB[inn,it_3+1]<agrid[it_3+1,1]
                    asimB[inn,it_3+1] = agrid[it_3+1,1]
                    csimB[inn,it_3] = xsimB[inn,it_3] - asimB[inn,it_3+1]
                end
            elseif it_3==Twork
                if asimB[inn,it_3+1]<agridret[1,1,round(Int,psimIB[inn])]
                    asimB[inn,it_3+1] = agridret[1,1,round(Int,psimIB[inn])]
                    csimB[inn,it_3] = xsimB[inn,it_3] - asimB[inn,it_3+1]
                end
            end
            #end if receiving shock
            if isnan(csim[inn,it_3])>0
                display("bad consumption")
            end

            it_3=it_3+1 # Go to next period
        end
    end
end

function Retirement_simulation(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)
    
    Threads.@threads for inn in 1:nsim
        for it in 1:Tret
            nnsim[inn, Twork+it] = inn;
            agesim[inn, Twork+it]= Twork+it;
            ysim[inn, Twork+it] = pgrid[it,round(Int,psimI[inn])];
            ypresim[inn, Twork+it] = ppregrid[it,round(Int, psimI[inn])]
            if it<Tret
                if asim[inn,Twork+it]<0
                    trsim[inn,Twork+it]    = max( cfloor - (Rnetdebt*asim[inn,Twork+it]+ ysim[inn,Twork+it] - agridret[it+1,1,round(Int, psimI[inn])]) ,0.0)
                else
                    trsim[inn,Twork+it]    = max( cfloor - (Rnetsave*asim[inn,Twork+it]+ ysim[inn,Twork+it] - agridret[it+1,1,round(Int, psimI[inn])]) ,0.0)
                end
                csim[inn,Twork+it] = LinInterp1(agridret[it,:,round(Int,psimI[inn])],conret[it,:,round(Int,psimI[inn])],asim[inn,Twork+it])
                csimT[inn,Twork+it] = (R-1).*(asim[inn,Twork+it])+ysim[inn,Twork+it];
                if asim[inn,Twork+it]<0
                    xsim[inn,Twork+it] = Rnetdebt*asim[inn,Twork+it]+ ysim[inn,Twork+it]  + trsim[inn,Twork+it]
                else
                    xsim[inn,Twork+it] = Rnetsave*asim[inn,Twork+it]+ ysim[inn,Twork+it]  + trsim[inn,Twork+it]
                end
                asim[inn,Twork+it+1] = (xsim[inn,Twork+it]  - csim[inn,Twork+it])/annprem[it]
            end
            if it<Tret
                if asim[inn,Twork+it+1]<agridret[it+1,1,round(Int,psimI[inn])]
                    asim[inn,Twork+it+1] = agridret[it+1,1,round(Int,psimI[inn])]
                    csim[inn,Twork+it] = xsim[inn,Twork+it] - asim[inn,Twork+it+1]*annprem[it]
                end
            end
        end
    end
end

function Result_simulation(;nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)

    alsim = zeros(nsim, Ttot)
    alsim = algrid[alsimI[:]]*ones(Ttot)'
    
    # Mean values for all simulations
    csim_mean = zeros(Ttot);
    csimB_mean = zeros(Ttot);
    csimT_mean = zeros(Ttot);
    asim_mean = zeros(Ttot);
    ysim_mean = zeros(Ttot);
    ysimB_mean = zeros(Ttot);
    
    for age in 1:(Ttot-1)
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

    # Mean values for all simulations
    result = Dict(
        :csim       => csim,
        :csim_mean  => csim_mean,
        
        :csimB      => csimB,
        :csimB_mean => csimB_mean,
        
        :csimT      => csimT,
        :csimT_mean => csimT_mean,

        :asim       => asim,
        :asim_mean  => asim_mean,
        
        :ysim       => ysim,
        :ysim_mean  => ysim_mean,

        :ysimB      => ysimB,
        :ysimB_mean => ysimB_mean,
        
        :zsim      => zsim,
        :zsim_mean => zsim
    )

    return result
end