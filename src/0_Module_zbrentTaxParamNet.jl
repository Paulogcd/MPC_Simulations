"""
Missing documentation
"""
function zbrentTaxParamNet(x1, x2, tol, func)

    e = zeros(1)
    d = zeros(1)

    itmax=200
    epsilon=eps(x1)

    a=x1
    b=x2
    fa=func(a)
    fb=func(b)

    if sign.(fa[1]*fb[1])>0
        display("root must be bracketed for zbrent")
    end

    c=b
    fc=fb
    iter=1
    while iter<=itmax

        if sign.(fb[1]*fc[1])>0
            c=a
            fc=fa
            d=b-a
            e=d
            iter=iter+1
        end

        if abs(fc[1])<abs(fb[1])
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        end

        tol1=2.0*epsilon*abs(b)+0.5*tol
        xm=0.5*(c-b)

        if (abs(xm)<=tol1)||(fb[1]==0)
            zbrent=b
            fb, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost = func(b)

            return [zbrent, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost]; #Just to end statement
        end

        if (abs(e)>=tol1)&&(abs(fa[1])>abs(fb[1]))

            s=fb[1]/fa[1]

            if a==c
                p=2.0*xm*s
                q=1.0-s
            else
                q=fa[1]/fc[1]
                r=fb[1]/fc[1]
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
                q=(q-1.0)*(r-1.0)*(s-1.0)
            end

            if p>0
                q=-q
            end

            p=abs(p)

            if (2.0*p < min(3.0*xm*q-abs(tol1*q), abs(e*q)))
                e=d
                d=p/q
            else
                d=xm
                e=d
            end

        else
            d=xm
            e=d
        end

        a=b
        fa=fb

        if abs(d)>tol1
            b=b+d
        else
            b=b+tol1*sign.(xm)
        end

        #d if abs(d)>tol1 and sign(tol1,xm) if abs(d)<=tol1
        fb, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost = func(b)
        iter=iter+1
    end

    display("zbrentTaxParam: exceeded maximum iterations")
    zbrent=b

    [zbrent, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost]
end

export zbrentTaxParamNet