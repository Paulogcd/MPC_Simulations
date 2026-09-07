"""
Missing documentation
"""
function rtsecSSParam(lx1, lx2, facc, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement)

    maxit=30
    x1 = lx1
    x2 = lx2

    fl=FnSSParam(lx1, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement)[1]
    if abs(fl)<facc
        disp("early return 1()")
        return [rtsec, pgrid, ppregrid, pind]
    end

    f=FnSSParam(lx2, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement)[1]
    if abs(f)<facc
        disp("early return 2()")
        return [rtsec, pgrid, ppregrid, pind]
    end

    if abs(fl) < abs(f)
        rtsec=x1
        xl=x2
        ltemp = fl
        fl = f
        f = ltemp
    else
        ()
        xl=x1
        rtsec=x2
    end

    j=1
    while j<=maxit
        dx=(xl-rtsec)*f/(f-fl)
        xl=rtsec
        fl=f
        rtsec=rtsec+dx
        f, pgrid, ppregrid, pind=FnSSParam(rtsec, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement)
        #f=func[rtsec]
        if abs(fl)<facc
            display("early return 3 ()")
            display(j)
            return [rtsec, pgrid, ppregrid, pind]
        end
        j=j+1
    end

    [rtsec, pgrid, ppregrid, pind]

end