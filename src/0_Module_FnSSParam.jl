
"""
Missing documentation
"""
function FnSSParam(lsspar, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement)

    # INTEREST RATE

    R = 1.03;                           # annual gross interest rate

    # PREFERENCE PARAMETERS 

    gam = 1.0;
    bet = 1/R;
    qprefb = 120000.0; # 10000000.0, bliss point
    qprefa = 1.0; # curvature in pref 

    # BORROWING LIMIT: SET TO VERY LARGE NEGATIVE NO. FOR NBL

    borrowlim = -100000000;             # 0.0

    # GOVERNMENT PARAMETERS 

    # gouveia strauss 
    stax = 2.0e-4;                      # guess: is chosen optimally
    ptax = 0.768;
    btax = 0.258;

    # other taxes and benefits
    cfloor = -100000000000.0;
    pentax = 0.0;    # payroll tax   
    rtax = 0.0; # tax on interest income
    otax = 0.0; # tax on old age pensions, check option in Parameters
    pencapfrac = 2.2; # cap on [pre-tax] earnings that contribute to pension index, as fraction of av pre-tax earnigs

    Rnet = 1.0 + (1.0-rtax)*(R-1);    # annual after tax interest rate


    g0=2.746
    g1=0.0624
    g2=-0.00167


    a=-2.495
    b=-0.1037
    c=-5.051
    d=-0.1087

    pind = zeros(ngpm, ngpz, ngpe, ngpal)
    ppregrid = zeros(Tret, ngpal*ngpe*ngpz*ngpm)
    pgrid = zeros(Tret, ngpal*ngpe*ngpz*ngpm)


    if Display==1
        display(" SS benefit parameter: ")
        display(lsspar)
    end

    if BendPointsPostTax==1
        lavearns = sum(avearnspost)/real(Twork)
    else
        lavearns = sum(avearnspre)/real(Twork)
    end

    if ScaleBendPoints==1
        ;
        ssbend1 = 0.18*lsspar*lavearns
        ssbend2 = 1.1*lsspar*lavearns
    else
        ssbend1 = 0.18*lavearns
        ssbend2 = 1.1*lavearns
    end


    ltotpen = 0
    ip=1
    im=1
    while im<=ngpm
        iz=1
        while iz<=ngpz
            ie=1
            while ie<=ngpe
                ia=1
                while ia<=ngpal
                    pind[im, iz, ie, ia] = ip
                    it=1
                    while it<=Tret
                        if UseFinalZPension==1
                            g=g0+g1*Twork+g2*Twork*Twork
                            xi=a+b*Twork+c*zgrid[Twork, iz]+d*it*zgrid[Twork, iz]
                            #I consider last income OR the last income that would be earned absent unemployment (else last income to compute benefist is zero)
                            lmearns = exp(g + algrid[ia] + egrid[Twork, ie] + zgrid[Twork, iz])
                        else
                            lmearns = mgrid[Twork, im]
                        end
                        if lmearns <= ssbend1
                            ppregrid[it, ip] = 0.9*lmearns
                        elseif lmearns <= ssbend2
                            ppregrid[it, ip] = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
                        else
                            ppregrid[it, ip] = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
                        end
                        if ScaleBendPoints==0
                            ppregrid[it, ip] = ppregrid[it, ip]*lsspar
                        end
                        #!pgrid[it,ip] = ppregrid[it,ip]*(1.0-otax)
                        FnTax = btax*(ppregrid[it, ip]*0.85 - ((ppregrid[it, ip]*0.85)^(-ptax) + stax)^(-1.0/ptax)) + pentax*ppregrid[it, ip]*0.85
                        pgrid[it, ip] = max(ppregrid[it, ip] - FnTax, 1)
                        it=it+1
                    end
                    ip = ip+1
                    ia=ia+1
                end
                ie=ie+1
            end
            iz=iz+1
        end
        im=im+1
    end

    if MatchAggPensionBen==1
        if UseFinalZPension==1
            ltotpen = 0
            ia=1
            while ia<=ngpal
                iz=1
                while iz<=ngpz
                    ie=1
                    while ie<=ngpe
                        ip = pind[1, iz, ie, ia]
                        it=1
                        while it<=Tret
                            ltotpen = ltotpen + ppregrid[it, ip]*zdist[Twork, iz]*edist[Twork, ie]*adist[ia]*popsize[Twork+it]
                            it=it+1
                        end
                        ip=ip+1
                        ie=ie+1
                    end
                    iz=iz+1
                end
                ia=ia+1
            end
        else
            inm=1
            while inm<=nsim
                lmloc, lp=FindLinProb1(mgrid[Twork-1, :], yavsim[inm, Twork-1])
                im=random(2, lp, 1)
                lpsim[inm] = ppregrid[1, pind[lmloc[im], zsimI[inm, Twork], esimI[inm, Twork], 1]]
                inm=inm+1
            end
            ltotpen = 0
            it = 1
            while it<=Tret
                ltotpen = ltotpen + popsize[Twork+it]*sum(lpsim)/real(nsim)
                it=it+1
            end
        end
        if UseNetIncToScaleSS == 1
            FnSSParam = ltotpen/totlabincpost - targetSSBenToLabinc
            if Display==1
                display(" SS Ben / Post-tax labor income: ")
                display((ltotpen/totlabincpost)*100)

            end
        else
            FnSSParam = ltotpen/totlabincpre - targetSSBenToLabinc
            if Display==1
                display(" SS Ben / Pre-tax labor income: ")
                display((ltotpen/totlabincpre)*100)
            end
        end
    else

        lmearns = sum(min.(avearnspre[:], pencap))/real(Twork)
        if lmearns <= ssbend1
            lmpen = 0.9*lmearns
        elseif lmearns <= ssbend2
            lmpen = 0.9*ssbend1 + 0.32*(lmearns - ssbend1)
        else
            lmpen = 0.9*ssbend1 + 0.32*(ssbend2-ssbend1) + 0.15*(lmearns - ssbend2)
        end

        if ScaleBendPoints==0
            lmpen = lmpen*lsspar
        end
        FnSSParam = lmpen/(sum(avearnspre)/real(Twork)) - targetSSAvReplacement
        if Display==1
            display(" SS Ben for Mean AIME/ Av Pre-tax Earns: ")
            display((lmpen/(SUM[avearnspre]/real(Twork)))*100)
        end
    end

    [FnSSParam, pgrid, ppregrid, pind]

end