"""
Missing documentation
"""
function golden(ax, bx, cx, tol, func)


    r=0.61803399
    c=1.0-r

    x0=ax;  #!At any given time we will keep trace of 4 points: 
    x3=cx;  #!X0,X1,X2,X3. 

    if abs(cx-bx)>abs(bx-ax)
        x1=bx
        x2=bx+c*(cx-bx)
    else
        x2=bx
        x1=bx-c*(bx-ax)
    end

    f1=func(x1)
    f2=func(x2)

    while abs(x3-x0)>tol*(abs(x1)+abs(x2))
        if f2<f1
            x0=x1
            x1=x2
            x2=r*x1+c*x3
            f1=f2
            f2=func(x2)
        else
            ()
            x3=x2
            x2=x1
            x1=r*x2+c*x0
            f2=f1
            f1=func(x1)
        end
    end

    if f1<f2
        golden=f1
        xmin=x1
    else
        golden=f2
        xmin=x2
    end
    xmin
end
