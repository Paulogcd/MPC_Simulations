export AgesResolution

"""
Solve for each age.
Equivalent of the former Module_DecisionsSImulate30Nucbar_Parallel.jl script.
"""
function AgesResolution(;Ttot, xgridret, ngpp, qprefa, qprefb, cbar, gam, annprem, pgrid)

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #! Start with last period; eat everything  !!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if Display == 1
        display("Solving for decision rules at age ")
        display(Ttot)
    end;

    conret = zeros(size(xgridret));
    mucret = zeros(size(xgridret));

    ip = 1;
    while ip <= ngpp
        conret[Tret, :, ip] = xgridret[Tret, :, ip]
        if QuadraticPref == 1
            mucret[Tret, :, ip] = qprefa .* (qprefb .- conret[Tret, :, ip])
        else
            mucret[Tret, :, ip] = (conret[Tret, :, ip] .- cbar) .^ (- gam[Tret]);
        end
        ip = ip+1
    end;

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #!! Retired Agents: Pension value as state variable !!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    mucret1 = zeros(size(mucret));
    conret1 = zeros(size(mucret));
    assret1 = zeros(size(mucret));

    size_a, size_b, size_c = size(mucret);
    assret = zeros((size_a-1, size_b, size_c));

    it = Tret - 1;

    while it >= 1
        if Display == 1
            display("Solving for decision rules at age ")
            display(Twork+it)
        end

        # PARALLEL

        @threads for ip in 1:ngpp

            # solve on tt+1 grid()
            if minimum(agridret[it+1, :, ip])>=0
                BLbind = 0
                mucret1[it, :, ip] = bet[Twork+it]*(surprob[it]/annprem[it])*Rnetsave*mucret[it+1, :, ip]
                if QuadraticPref==1
                    conret1[it, :, ip] = qprefb .- mucret1[it, :, ip] ./ qprefa
                else
                    conret1[it, :, ip] = mucret1[it, :, ip] .^ (-1.0/gam[it]) .+ cbar;
                end
                assret1[it, :, ip] = (conret1[it, :, ip] + agridret[it+1, :, ip] * annprem[it] .- pgrid[it, ip] - tranret[it, :, ip])/Rnetsave
            else
                A = agridret[it+1, :, ip]
                B, BLbind = findmax(A[A .< 0])
                BLbind = round(Int, BLbind)
                mucret1[it, 1:BLbind, ip] = bet[Twork+it]*(surprob[it]/annprem[it])*Rnetdebt*mucret[it+1, 1:BLbind, ip]
                mucret1[it, (BLbind+1):ngpa, ip] = bet[Twork+it] * (surprob[it]/annprem[it]) * Rnetsave * mucret[it+1, (BLbind+1):ngpa, ip]
                if QuadraticPref == 1
                    conret1[it, :, ip] = qprefb .- mucret1[it, :, ip] ./ qprefa
                else
                    conret1[it, :, ip] = mucret1[it, :, ip] .^ (-1.0/gam[it]) .+ cbar;
                end
                assret1[it, 1:BLbind, ip] = (conret1[it, 1:BLbind, ip]+agridret[it+1, 1:BLbind, ip]*annprem[it] .- pgrid[it, ip]-tranret[it, 1:BLbind, ip])/Rnetdebt
                assret1[it, (BLbind+1):ngpa, ip] = (conret1[it, (BLbind+1):ngpa, ip]+agridret[it+1, (BLbind+1):ngpa, ip]*annprem[it] .- pgrid[it, ip] - tranret[it, (BLbind+1):ngpa, ip])/Rnetsave
            end

            # deal with borrowing limits
            if minimum(agridret[it, :, ip])>=assret1[it, 1, ip]
                BLbind = 0
            else
                A = agridret[it, :, ip]
                B, BLbind = findmax(A[A .< assret1[it, 1, ip]])
                BLbind = round(Int, BLbind)
                assret[it, 1:round(Int, BLbind), ip] .= agridret[it+1, 1, ip]
                conret[it, 1:round(Int, BLbind), ip] = xgridret[it, 1:round(Int, BLbind), ip] - assret[it, 1:round(Int, BLbind), ip]*annprem[it]
                if QuadraticPref==1
                    mucret[it, 1:round(Int, BLbind), ip] = qprefa .* (qprefb .- conret[it, 1:round(Int, BLbind), ip])
                else
                    mucret[it, 1:round(Int, BLbind), ip] = (conret[it, 1:round(Int, BLbind), ip] .- cbar) .^ (-gam[it])
                end
            end

            # interpolate muc1 as fun of ass1 to get muc where BL does not bind

            conret[it, (BLbind+1):ngpa, ip] = LinInterp(ngpa, assret1[it, :, ip], conret1[it, :, ip], ngpa-BLbind, agridret[it, (BLbind+1):ngpa, ip])
            if QuadraticPref==1
                mucret[it, (round(Int, BLbind)+1):ngpa, ip] = qprefa .* (qprefb .- conret[it, (round(Int, BLbind)+1):ngpa, ip])
            else
                mucret[it, (round(Int, BLbind)+1):ngpa, ip] = (conret[it, (round(Int, BLbind)+1):ngpa, ip] .- cbar) .^ (-gam[it])
            end
            assret[it, (round(Int, BLbind)+1):ngpa, ip] = (xgridret[it, (round(Int, BLbind)+1):ngpa, ip] - conret[it, (round(Int, BLbind)+1):ngpa, ip])/annprem[it]

            if any(isnan, conret[it, :, ip])>0
                display("nan encountered")
            end

        end
        @info(string("Age: ", Twork+it, ". Time retirement period: "));
        it=it-1
    end;


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
    #!! Working Agents: Rules depend on shocks  !!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    muc1 = zeros(Twork, ngpa, ngpm, ngpz, ngpe, ngpal, 2);
    con1 = zeros(Twork, ngpa, ngpm, ngpz, ngpe, ngpal, 2);
    ass1 = zeros(Twork, ngpa, ngpm, ngpz, ngpe, ngpal, 2);
    muc = zeros(Twork, ngpa, ngpm, ngpz, ngpe, ngpal, 2);
    con = zeros(Twork, ngpa, ngpm, ngpz, ngpe, ngpal, 2);
    ass = zeros(Twork, ngpa, ngpm, ngpz, ngpe, ngpal, 2);

    it = Twork;
    while it >= 1
        if Display == 1
            @info(string("Solving for decision rules at age ", it, "..."))
        end
        inu = 1
        while inu <= 2
            ial = 1
            while ial <= ngpal
                @threads for iz in 1:ngpz
                    @threads for ie in 1:ngpe
                        im=1
                        while im <= ngpm
                            ip = round(Int, pind[im, iz, ie, ial])
                            if it == Twork
                                if minimum(agridret[1, :, ip]) >= 0
                                    ALbind = 0
                                    muc1[it, :, im, iz, ie, ial, inu] = bet[it]*Rnetsave*mucret[1, :, round(Int, pind[im, iz, ie, ial])];# Euler equation
                                    if QuadraticPref==1
                                        con1[it, :, im, iz, ie, ial, inu] = qprefb .- muc1[it, :, im, iz, ie, ial, inu] ./ qprefa
                                    else
                                        con1[it, :, im, iz, ie, ial, inu] = muc1[it, :, im, iz, ie, ial, inu] .^ (-1.0/gam[it]) .+ cbar;
                                    end
                                    ass1[it, :, im, iz, ie, ial, inu] = (con1[it, :, im, iz, ie, ial, inu] + agridret[1, :, ip] .- ygrid[it, iz, ie, ial, inu] - tran[it, :, iz, ie, ial, inu])/Rnetsave;
                                    if minimum(agrid[it, :]) >= ass1[it, 1, im, iz, ie, ial, inu]#!BL does not bind anywhere
                                        BLbind=0
                                    else
                                        A=agrid[it, :]
                                        B, BLbind=findmax(A[A .< ass1[it, 1, im, iz, ie, ial, inu]])
                                        ass[it, 1:BLbind, im, iz, ie, ial, inu] .= agridret[1, 1, ip]
                                        con[it, 1:BLbind, im, iz, ie, ial, inu] = xgrid[it, 1:BLbind, iz, ie, ial, inu]-ass[it, 1:BLbind, im, iz, ie, ial, inu]
                                        if QuadraticPref==1
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, 1:BLbind, im, iz, ie, ial, inu])
                                        else
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = (con[it, 1:BLbind, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                        end
                                    end
                                    con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = LinInterp(ngpa, ass1[it, :, im, iz, ie, ial, inu], con1[it, :, im, iz, ie, ial, inu], ngpa-BLbind, agrid[it, (BLbind+1):ngpa])
                                    if QuadraticPref==1
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu])
                                    else
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = (con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                    end
                                    ass[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = xgrid[it, (BLbind+1):ngpa, iz, ie, ial, inu] - con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu];#!budget constraint
                                else
                                    A = agridret[1, :, ip]
                                    A, ALbind = findmax(A[A .< 0])
                                    ALbind = round(Int, ALbind)
                                    muc1[it, 1:ALbind, im, iz, ie, ial, inu] = bet[it]*Rnetdebt*mucret[1, 1:ALbind, round(Int, pind[im, iz, ie, ial])];
                                    muc1[it, (ALbind+1):ngpa, im, iz, ie, ial, inu] = bet[it]*Rnetsave*mucret[1, (ALbind+1):ngpa, round(Int, pind[im, iz, ie, ial])];
                                    if QuadraticPref==1
                                        con1[it, :, im, iz, ie, ial, inu] = qprefb .- muc1[it, :, im, iz, ie, ial, inu] ./ qprefa
                                    else
                                        con1[it, :, im, iz, ie, ial, inu] = muc1[it, :, im, iz, ie, ial, inu] .^ (-1.0/gam[it]) .+ cbar;
                                    end
                                    ass1[it, 1:ALbind, im, iz, ie, ial, inu] = (con1[it, 1:ALbind, im, iz, ie, ial, inu] + agridret[1, 1:ALbind, ip] .- ygrid[it, iz, ie, ial, inu] - tran[it, 1:ALbind, iz, ie, ial, inu])/Rnetdebt;
                                    ass1[it, (ALbind+1):ngpa, im, iz, ie, ial, inu] = (con1[it, (ALbind+1):ngpa, im, iz, ie, ial, inu] + agridret[1, (ALbind+1):ngpa, ip] .- ygrid[it, iz, ie, ial, inu] - tran[it, (ALbind+1):ngpa, iz, ie, ial, inu])/Rnetsave;
                                    if minimum(agrid[it, :]) >= ass1[it, 1, im, iz, ie, ial, inu]#!BL does not bind anywhere
                                        BLbind=0
                                    else
                                        A = agrid[it, :]
                                        B, BLbind = findmax(A[A .< ass1[it, 1, im, iz, ie, ial, inu]])
                                        ass[it, 1:BLbind, im, iz, ie, ial, inu] .= agridret[1, 1, ip]
                                        con[it, 1:BLbind, im, iz, ie, ial, inu] = xgrid[it, 1:BLbind, iz, ie, ial, inu]-ass[it, 1:BLbind, im, iz, ie, ial, inu]
                                        if QuadraticPref == 1
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, 1:BLbind, im, iz, ie, ial, inu])
                                        else
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = (con[it, 1:BLbind, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                        end
                                    end
                                    con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = LinInterp(ngpa, ass1[it, :, im, iz, ie, ial, inu], con1[it, :, im, iz, ie, ial, inu], ngpa-BLbind, agrid[it, (BLbind+1):ngpa])
                                    if QuadraticPref==1
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu])
                                    else
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = (con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                    end
                                    ass[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = xgrid[it, (BLbind+1):ngpa, iz, ie, ial, inu] - con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu];#!budget constraint
                                end
                            else
                                emuc = zeros(size(muc[it, :, im, iz, ie, ial, inu]))
                                iz2=1
                                while iz2<=ngpz
                                    ie2=1
                                    while ie2<=ngpe
                                        # find two probabilities and average over:
                                        lnextm = (it*mgrid[it, im] + min(ypregrid[it+1, iz2, ie2, ial, inu], pencap))/real(it+1)
                                        imnext, lpmnext=FindLinProb1(mgrid[it+1, :], lnextm)
                                        imnext = round.(Int, imnext)
                                        p=pzgrid[it, iz2]
                                        emuc[:] = emuc[:] + (lpmnext[1]*muc[it+1, :, imnext[1], iz2, ie2, ial, 1]*ztrans[it, iz, iz2]*edist[it+1, ie2] + lpmnext[2]*muc[it+1, :, imnext[2], iz2, ie2, ial, 1]*ztrans[it, iz, iz2]*edist[it+1, ie2])*(1-p) + (lpmnext[1]*muc[it+1, :, imnext[1], iz2, ie2, ial, 2]*ztrans[it, iz, iz2]*edist[it+1, ie2] + lpmnext[2]*muc[it+1, :, imnext[2], iz2, ie2, ial, 2]*ztrans[it, iz, iz2]*edist[it+1, ie2])*(p);
                                        ie2=ie2+1
                                    end
                                    iz2=iz2+1
                                end
                                if minimum(agrid[it+1, :])>=0 #t<Twork and R=Rsave
                                    ALbind=0
                                    muc1[it, :, im, iz, ie, ial, inu] = bet[it]*Rnetsave*emuc; # Euler equation
                                    if QuadraticPref==1
                                        con1[it, :, im, iz, ie, ial, inu] = qprefb .- muc1[it, :, im, iz, ie, ial, inu] ./ qprefa
                                    else
                                        con1[it, :, im, iz, ie, ial, inu] = muc1[it, :, im, iz, ie, ial, inu] .^ (-1.0/gam[it]) .+ cbar;
                                    end
                                    ass1[it, :, im, iz, ie, ial, inu] = (con1[it, :, im, iz, ie, ial, inu] + agrid[it+1, :] .- ygrid[it, iz, ie, ial, inu] - tran[it, :, iz, ie, ial, inu])/Rnetsave;
                                    if minimum(agrid[it, :]) >= ass1[it, 1, im, iz, ie, ial, inu]#!BL does not bind anywhere
                                        BLbind=0
                                    else
                                        A=agrid[it, :]
                                        B, BLbind=findmax(A[A .< ass1[it, 1, im, iz, ie, ial, inu]])
                                        if it==Twork
                                            ass[it, 1:BLbind, im, iz, ie, ial, inu] .= agridret[1, 1, ip]
                                        else
                                            ass[it, 1:BLbind, im, iz, ie, ial, inu] .= agrid[it+1, 1]
                                        end
                                        con[it, 1:BLbind, im, iz, ie, ial, inu] = xgrid[it, 1:BLbind, iz, ie, ial, inu]-ass[it, 1:BLbind, im, iz, ie, ial, inu]
                                        if QuadraticPref==1
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, 1:BLbind, im, iz, ie, ial, inu])
                                        else
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = (con[it, 1:BLbind, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                        end
                                    end
                                    con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = LinInterp(ngpa, ass1[it, :, im, iz, ie, ial, inu], con1[it, :, im, iz, ie, ial, inu], ngpa-BLbind, agrid[it, (BLbind+1):ngpa])
                                    if QuadraticPref==1
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu])
                                    else
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = (con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                    end
                                    ass[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = xgrid[it, (BLbind+1):ngpa, iz, ie, ial, inu] - con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu];#!budget constraint
                                else #t<Twork and R=Rsave or Rdebt
                                    A=agrid[it+1, :]
                                    B, ALbind=findmax(A[A .< 0])
                                    ALbind = round(Int, ALbind)
                                    muc1[it, 1:ALbind, im, iz, ie, ial, inu] = bet[it]*Rnetsave*emuc[1:ALbind];
                                    muc1[it, (ALbind+1):ngpa, im, iz, ie, ial, inu] = bet[it]*Rnetsave*emuc[(ALbind+1):ngpa];
                                    if QuadraticPref==1
                                        con1[it, :, im, iz, ie, ial, inu] = qprefb .- muc1[it, :, im, iz, ie, ial, inu] ./ qprefa
                                    else
                                        con1[it, :, im, iz, ie, ial, inu] = muc1[it, :, im, iz, ie, ial, inu] .^ (-1.0/gam[it]) .+ cbar;
                                    end
                                    ass1[it, 1:ALbind, im, iz, ie, ial, inu] = (con1[it, 1:ALbind, im, iz, ie, ial, inu] + agrid[it+1, 1:ALbind] .- ygrid[it, iz, ie, ial, inu] - tran[it, 1:ALbind, iz, ie, ial, inu])/Rnetdebt
                                    ass1[it, (ALbind+1):ngpa, im, iz, ie, ial, inu] = (con1[it, (ALbind+1):ngpa, im, iz, ie, ial, inu] + agrid[it+1, (ALbind+1):ngpa] .- ygrid[it, iz, ie, ial, inu] - tran[it, (ALbind+1):ngpa, iz, ie, ial, inu])/Rnetsave
                                    if minimum(agrid[it, :]) >= ass1[it, 1, im, iz, ie, ial, inu]#!BL does not bind anywhere
                                        BLbind=0
                                    else
                                        A = agrid[it, :]
                                        B, BLbind = findmax(A[A .< ass1[it, 1, im, iz, ie, ial, inu]])
                                        if it == Twork
                                            ass[it, 1:BLbind, im, iz, ie, ial, inu] .= agridret[1, 1, ip]
                                        else
                                            ass[it, 1:BLbind, im, iz, ie, ial, inu] .= agrid[it+1, 1]
                                        end
                                        con[it, 1:BLbind, im, iz, ie, ial, inu] = xgrid[it, 1:BLbind, iz, ie, ial, inu]-ass[it, 1:BLbind, im, iz, ie, ial, inu]
                                        if QuadraticPref==1
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, 1:BLbind, im, iz, ie, ial, inu])
                                        else
                                            muc[it, 1:BLbind, im, iz, ie, ial, inu] = (con[it, 1:BLbind, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                        end
                                    end
                                    con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = LinInterp(ngpa, ass1[it, :, im, iz, ie, ial, inu], con1[it, :, im, iz, ie, ial, inu], ngpa-BLbind, agrid[it, (BLbind+1):ngpa])
                                    if QuadraticPref==1
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = qprefa .* (qprefb .- con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu])
                                    else
                                        muc[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = (con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] .- cbar) .^ (-gam[it])
                                    end
                                    ass[it, (BLbind+1):ngpa, im, iz, ie, ial, inu] = xgrid[it, (BLbind+1):ngpa, iz, ie, ial, inu] - con[it, (BLbind+1):ngpa, im, iz, ie, ial, inu];#!budget constraint
                                end
                            end
                            im=im+1
                        end
                    end
                end
                ial = ial + 1
            end
            inu = inu + 1
        end;
        @info(string("Age: ", it, ". Time retirement period: "));
        it = it-1
    end;

    resultat = (; conret, mucret, conret1, mucret1, assret1, muc1, con1, ass1, muc, con, ass)
    return resultat;
end