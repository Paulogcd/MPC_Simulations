export CashInHands

"""
Replace the long in-book function and produces the Cash-in-hands grids.
"""
function CashInHands(;Twork, Tret, ngpa, ngpz, ngpe, ngpal, ngpp, ygrid, agrid, cfloor, Rnetsave, Rnetdebt, pind)

    tran = zeros(Twork, ngpa, ngpz, ngpe, ngpal, 2);
    xgrid = zeros(Twork, ngpa, ngpz, ngpe, ngpal, 2);
    tranret = zeros(Tret, ngpa, ngpp);
    xgridret = zeros(Tret, ngpa, ngpp);

    ia = 1;

    while ia <= ngpa
        ial = 1
        while ial <= ngpal
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
            ial = ial + 1
        end
        ia = ia + 1
    end;

    resultat = (;tran, xgrid, tranret, xgridret)

    return resultat
end