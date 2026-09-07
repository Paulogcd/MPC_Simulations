"""
Missing documentation
"""
function mnbrak(ax, bx, func)

    # SUBROUTINE MNBRAK[AX,BX,CX,FA,FB,FC,FUNC]
    # !Given a function FUNC[X], and given distinct initial points AX and
    # !BX, this routine searches in the downhill direction (defined by the
    # !function as evaluated at the initial points) and returns new points
    # !AX, BX, CX which bracket a minimum of the function(). Also returned
    # !are the function values at the three points, FA, FB and FC.

    gold=1.618034
    glimit=100
    tiny=1.0e-20

    fa=func(ax)
    fb=func(bx)
    if fb > fa
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
    end

    cx=bx+gold*(bx-ax)
    fc=func(cx)
    iter=1

    while (fb>=fc)&&(iter<=25)
        Proceed = true
        iter=iter+1
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0 * (sign.(q-r)) * max(abs(q-r), tiny))
        ulim=bx+glimit*(cx-bx)
        if ((bx-u)*(u-cx) > 0)&&(Proceed)
            fu=func(u)
            if fu < fc
                ax=bx
                fa=fb
                bx=u
                fb=fu
                Proceed=false
                break
            elseif fu > fb
                cx=u
                fc=fu
                Proceed=false
                break
            end
            u=cx+gold*(cx-bx)
            fu=func(u)
        elseif (cx-u)*(u-ulim) > 0
            fu=func(u)
            if fu < fc
                bx=cx
                cx=u
                u=cx+gold*(cx-bx)
                fb=fc
                fc=fu
                fu=func(u)
            end
        elseif (u-ulim)*(ulim-cx) >= 0
            u=ulim
            fu=func(u)
        else
            u=cx+gold*(cx-bx)
            fu=func(u)
        end
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
    end


    if (iter>=100)
        cx=ax
        fc=fa
    end
    axx=ax
    bxx=bx

    [axx, bxx, cx, fa, fb, fc]
end