export Grids_formation

"""
Composed and intermediate function creating the grids used later in the model resolution.
"""
function Grids_formation(;nsim, Twork, rng, ngpz, ngpe, zdist, edist, ngpal, aldist, Ttot, pzgrid, ypregrid, ngpm, Tret, kappa, algrid, egrid, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement, ztrans, Pensions, ngpp, Rnetdebt, annprem, ygrid, ngpa, borrowlim, pexpgrid, amax, cfloor, Rnetsave);

    result_A = Grid_formation_A(;nsim, Twork, rng, ngpz, ngpe, zdist, edist, ngpal, aldist, Ttot, pzgrid, ypregrid, pencap, ztrans);

    result_B = Grid_formation_B(;nsim, Twork, ngpm, result_A.yavsim, Tret, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement, Pensions, ngpp);

    result_C = Grid_formation_3_natural_borrowing_constraints(;Tret, ngpp, Rnetdebt, annprem, result_B.pgrid);

    # Working periods
    result_4 = Grid_formation_4_working_periods(;Twork, Tret, ngpal, result_C.nblret, Rnetdebt, result_B.pind, ygrid, ngpa, borrowlim, pexpgrid, amax, ngpp);

    # Cash in Hand
    result_5 = Grid_formation_5_cash_in_hands(;Twork, Tret, ngpa, ngpz, ngpe, ngpal, ngpp, ygrid, result_4.agrid, result_4.agridret, cfloor, Rnetsave, Rnetdebt, result_B.pind, result_B.pgrid);

    @info("Finished forming grids")

    result = (;
        result_A.ypresim,
        result_A.alsimI,
        result_A.zsimI,
        result_A.esimI,
        result_A.yavsim,

        result_B.pgrid,
        result_B.pind,
        result_B.mgrid,
        result_B.ppregrid,

        result_4.agrid,
        result_4.agridret,

        result_5.tran,
        result_5.xgrid,
        result_5.tranret,
        result_5.xgridret
    )

    return result

end

function Grid_formation_A(;nsim, Twork, rng, ngpz, ngpe, zdist, edist, ngpal, aldist, Ttot, pzgrid, ypregrid, pencap, ztrans);
    
    zsimI = zeros(Int8, nsim, Twork);
    esimI = zeros(Int8, nsim, Twork);
    alsimI = zeros(Int8, nsim);
    nusimI = zeros(Int8, nsim, Twork);

    # Initial period
    zsimI[:, 1] = sample(rng, collect(1:ngpz), weights(zdist[1, :]), nsim);
    esimI[:, 1] = sample(rng, collect(1:ngpe), weights(edist[1, :]), nsim);
    alsimI[:]   = sample(rng, collect(1:ngpal), weights(aldist[:]), nsim);

    ypresim = zeros(nsim, Ttot);
    yavsim = zeros(nsim, Twork);

    Threads.@threads for inn in 1:nsim

        nudist = Binomial(1, pzgrid[1, zsimI[inn, 1]])
        nu = rand(rng, nudist, 1)[1] # ν = 1 means unemployment; ν = 0 means employment
        
        #I compute only on employed people=#
        ypresim[inn, 1] = ypregrid[1, zsimI[inn, 1], esimI[inn, 1], alsimI[inn], 1]
        yavsim[inn, 1] = min(ypresim[inn, 1], pencap)
        
        for it in 2:Twork
            zsimI[inn, it] = sample(rng, collect(1:ngpz), weights(ztrans[it-1, zsimI[inn, it-1], :]), 1)[1, 1]
            esimI[inn, it] = sample(rng, collect(1:ngpe), weights(edist[it, :]), 1)[1, 1]
            nudist = Binomial(1, pzgrid[it, zsimI[inn, it]])
            nu = rand(rng, nudist, 1)[1]
            #I compute only on employed people=#
            ypresim[inn, it] = ypregrid[it, zsimI[inn, it], esimI[inn, it], alsimI[inn], 1]
            yavsim[inn, it] = ((it-1)*yavsim[inn, it-1] + min(ypresim[inn, it], pencap))/real(it)
        end

    end;
    result = (;zsimI, esimI, alsimI, nusimI, ypresim, yavsim)
    return result
end

function Grid_formation_B(;nsim, Twork, ngpm, yavsim, Tret, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement, Pensions, ngpp);
    # equally spaced between 5th perc and median(), and median and 95th perc

    mgrid = zeros(Twork, ngpm);

    Threads.@threads for it in 1:Twork

        lyavsim = yavsim[:, it]
        lyavsim = sort(lyavsim)
        lperc5  = lyavsim[round(Int, nsim*0.05)]
        lperc50 = lyavsim[round(Int, nsim*0.5)]
        lperc95 = lyavsim[round(Int, nsim*0.95)]

        mgrid[it, 1] = lperc5
        mgrid[it, round(Int, (1+ngpm)/2)] = lperc50
        mgrid[it, ngpm] = lperc95

        lwidth_1 = (mgrid[it, round(Int, (1+ngpm)/2)] - mgrid[it, 1])/real((ngpm-1)/2)
        lwidth2 = (mgrid[it, ngpm] - mgrid[it, round(Int, (1+ngpm)/2)])/real((ngpm-1)/2)

        for im in Int.(2:((ngpm-1)/2))
            mgrid[it, im] = mgrid[it, im-1] + lwidth_1
        end

        for im in (round(Int, (ngpm+1)/2 + 1)):(ngpm-1)
            mgrid[it, im] = mgrid[it, im-1] + lwidth2
        end

    end;

    # Pensions

    if Pensions==1
        rtsec, pgrid, ppregrid, pind = rtsecSSParam_lowr(0.1, 1.5, 1.0e-3, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement);
    else
        rtsec    = 1
        pgrid    = zeros(Tret, ngpp)
        ppregrid = zeros(Tret, ngpp)
        pind     = ones(ngpm, ngpz, ngpe, ngpal)
    end

    result = (;rtsec, pgrid, ppregrid, pind, mgrid)

    return result
end

function Grid_formation_3_natural_borrowing_constraints(;Tret, ngpp, Rnetdebt, annprem, pgrid)
    
    # Natural borrowing limits
    
    # Last period
    nblret = zeros(Tret, ngpp);

    Threads.@threads for ip in 1:ngpp
        nblret[Tret, ip] = -pgrid[Tret, ip]/Rnetdebt + 0.0001;
    end;

    # Retired
    for it in Int.((Tret-1):-1:1)
        for ip in 1:ngpp
            nblret[it, ip] = annprem[it]*nblret[it+1, ip]/Rnetdebt-pgrid[it, ip]/Rnetdebt + 0.0001;
        end
    end;
    
    result = (;nblret)
    return result

end

function Grid_formation_4_working_periods(;Twork, Tret, ngpal, nblret, Rnetdebt, pind, ygrid, ngpa, borrowlim, pexpgrid, amax, ngpp)
    
    # Last Working Period
    nbl = zeros(Twork, ngpal);
    Threads.@threads for ial in 1:ngpal
        nbl[Twork, ial] = nblret[1, round(Int, pind[1, 1, 1, ial])]/Rnetdebt - ygrid[Twork, 1, 1, ial, 2]/Rnetdebt + 0.0001;
    end

    # Working

    for it in (Twork-1):-1:1
        for ial=1:ngpal
            nbl[it, ial] = nbl[it+1, ial]/Rnetdebt - ygrid[it, 1, 1, ial, 2]/Rnetdebt + 0.0001;
        end
    end;

    # Assets: 0 to amax - exponentially spaced grid() - pexpgrid is parametr

    agrid = zeros(Twork, ngpa);

    Threads.@threads for it in 1:Twork
        agrid[it, 1] = max(borrowlim, nbl[it, 1])
        for ia in 2:ngpa
            agrid[it, ia] = agrid[it, ia-1] + exp(pexpgrid*(ia-1))
        end
        ltemp1 = (amax-agrid[it, 1])/(agrid[it, ngpa]-agrid[it, 1])
        agrid[it, :] = agrid[it, 1] .+ ltemp1 .* (agrid[it, :] .- agrid[it, 1])
    end;

    agridret = zeros(Tret, ngpa, ngpp);

    Threads.@threads for it in 1:Tret
        for ip in 1:ngpp
            agridret[it, 1, ip] = max(borrowlim, nblret[it, ip])
            for ia in 2:ngpa
                agridret[it, ia, ip] = agridret[it, ia-1, ip] + exp(pexpgrid*(ia-1))
            end
            ltemp1 = (amax-agridret[it, 1, ip])/(agridret[it, ngpa, ip]-agridret[it, 1, ip])
            agridret[it, :, ip] = agridret[it, 1, ip] .+ ltemp1 .* (agridret[it, :, ip] .- agridret[it, 1, ip])
        end
    end;
    
    result = (;agrid, agridret)
    return result
end

function Grid_formation_5_cash_in_hands(;Twork, Tret, ngpa, ngpz, ngpe, ngpal, ngpp, ygrid, agrid, agridret, cfloor, Rnetsave, Rnetdebt, pind, pgrid)

    tran = zeros(Twork, ngpa, ngpz, ngpe, ngpal, 2);
    xgrid = zeros(Twork, ngpa, ngpz, ngpe, ngpal, 2);
    tranret = zeros(Tret, ngpa, ngpp);
    xgridret = zeros(Tret, ngpa, ngpp);

    Threads.@threads for ia in 1:ngpa
        for ial in 1:ngpal
            iz = 1
            while iz <= ngpz
                ie = 1
                while ie <= ngpe
                    it = 1
                    if agrid[it, ia] >= 0
                        tran[it, ia, iz, ie, ial, 1] = max(cfloor - (Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] - agrid[it+1, 1]), 0.0)
                        xgrid[it, ia, iz, ie, ial, 1] = Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] + tran[it, ia, iz, ie, ial, 1]
                        tran[it, ia, iz, ie, ial, 2] = max(cfloor - (Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] - agrid[it+1, 1]), 0.0)
                        xgrid[it, ia, iz, ie, ial, 2] = Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] + tran[it, ia, iz, ie, ial, 2]
                    else
                        tran[it, ia, iz, ie, ial, 1] = max(cfloor - (Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] - agrid[it+1, 1]), 0.0)
                        xgrid[it, ia, iz, ie, ial, 1] = Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] + tran[it, ia, iz, ie, ial, 1]
                        tran[it, ia, iz, ie, ial, 2] = max(cfloor - (Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] - agrid[it+1, 1]), 0.0)
                        xgrid[it, ia, iz, ie, ial, 2] = Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] + tran[it, ia, iz, ie, ial, 2]
                    end
                    it=2
                    while it<=Twork-1
                        if agrid[it, ia] >=0
                            tran[it, ia, iz, ie, ial, 1] = max(cfloor - (Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] - agrid[it+1, 1]), 0.0)
                            xgrid[it, ia, iz, ie, ial, 1] = Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] + tran[it, ia, iz, ie, ial, 1]
                            tran[it, ia, iz, ie, ial, 2] = max(cfloor - (Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] - agrid[it+1, 1]), 0.0)
                            xgrid[it, ia, iz, ie, ial, 2] = Rnetsave * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] + tran[it, ia, iz, ie, ial, 2]
                        else
                            tran[it, ia, iz, ie, ial, 1] = max(cfloor - (Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] - agrid[it+1, 1]), 0.0)
                            xgrid[it, ia, iz, ie, ial, 1] = Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 1] + tran[it, ia, iz, ie, ial, 1]
                            tran[it, ia, iz, ie, ial, 2] = max(cfloor - (Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] - agrid[it+1, 1]), 0.0)
                            xgrid[it, ia, iz, ie, ial, 2] = Rnetdebt * agrid[it, ia] + ygrid[it, iz, ie, ial, 2] + tran[it, ia, iz, ie, ial, 2]
                        end
                        it=it+1
                    end
                    if agrid[Twork, ia] >=0
                        tran[Twork, ia, iz, ie, ial, 1] = max(cfloor - (Rnetsave * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 1] - agridret[1, 1, round(Int, pind[1, iz, ie, ial])]), 0.0)
                        xgrid[Twork, ia, iz, ie, ial, 1] = Rnetsave * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 1] + tran[Twork, ia, iz, ie, ial, 1]
                        tran[Twork, ia, iz, ie, ial, 2] = max(cfloor - (Rnetsave * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 2] - agridret[1, 1, round(Int, pind[1, iz, ie, ial])]), 0.0)
                        xgrid[Twork, ia, iz, ie, ial, 2] = Rnetsave * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 2] + tran[Twork, ia, iz, ie, ial, 2]
                    else
                        tran[Twork, ia, iz, ie, ial, 1] = max(cfloor - (Rnetdebt * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 1] - agridret[1, 1, round(Int, pind[1, iz, ie, ial])]), 0.0)
                        xgrid[Twork, ia, iz, ie, ial, 1] = Rnetdebt * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 1] + tran[Twork, ia, iz, ie, ial, 1]
                        tran[Twork, ia, iz, ie, ial, 2] = max(cfloor - (Rnetdebt * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 2] - agridret[1, 1, round(Int, pind[1, iz, ie, ial])]), 0.0)
                        xgrid[Twork, ia, iz, ie, ial, 2] = Rnetdebt * agrid[Twork, ia] + ygrid[Twork, iz, ie, ial, 2] + tran[Twork, ia, iz, ie, ial, 2]
                    end
                    ie = ie+1
                end
                iz = iz+1
            end
            ip = 1
            while ip <= ngpp
                it = 1
                while it <= Tret-1
                    if agridret[it, ia, ip]>0
                        tranret[it, ia, ip] = max(cfloor - (Rnetsave * agridret[it, ia, ip] + pgrid[it, ip] - agridret[it+1, 1, ip]), 0.0)
                        xgridret[it, ia, ip] = Rnetsave * agridret[it, ia, ip] + pgrid[it, ip] + tranret[it, ia, ip]
                    else
                        tranret[it, ia, ip] = max(cfloor - (Rnetdebt * agridret[it, ia, ip] + pgrid[it, ip] - agridret[it+1, 1, ip]), 0.0)
                        xgridret[it, ia, ip] = Rnetdebt * agridret[it, ia, ip] + pgrid[it, ip] + tranret[it, ia, ip]
                    end
                    it = it+1
                end
                if agridret[Tret, ia, ip]>0
                    tranret[Tret, ia, ip] = max(cfloor - (Rnetsave * agridret[Tret, ia, ip] + pgrid[Tret, ip]), 0.0)
                    xgridret[Tret, ia, ip] = Rnetsave * agridret[Tret, ia, ip] + pgrid[Tret, ip] + tranret[Tret, ia, ip]
                else
                    tranret[Tret, ia, ip] = max(cfloor - (Rnetdebt * agridret[Tret, ia, ip] + pgrid[Tret, ip]), 0.0)
                    xgridret[Tret, ia, ip] = Rnetdebt * agridret[Tret, ia, ip] + pgrid[Tret, ip] + tranret[Tret, ia, ip]
                end
                ip = ip + 1
            end
        end;
    end;

    result = (;
        tran,
        xgrid,
        tranret,
        xgridret
    )
    return result
end
