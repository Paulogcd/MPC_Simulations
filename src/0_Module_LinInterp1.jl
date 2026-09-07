"""
Missing documentation
"""
function LinInterp1(x, y, xi)
    # !this does linear interpolation of [x,y] at points xi
    # !requires x to be sorted in ascending order
    # !extrapolates out of range

    n = length(x)

    if minimum(x)>=xi
        v, LocL=findmax(x)
    else
        v, LocL=findmax(x[xi .> x])
    end
    #    [v,LocL] =max(x[xi>x])
    if xi<=x[1]
        LocL = 1
    end
    if LocL>=n
        LocL=n-1
    end

    xL = x[LocL]
    xH = x[LocL+1]
    yL = y[LocL]
    yH = y[LocL+1]

    yi = yL + (xi-xL)*((yH-yL)/(xH-xL))

    yi

end