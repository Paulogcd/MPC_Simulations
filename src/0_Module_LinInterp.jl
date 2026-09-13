"""
    Returns the linear interpolation of [x, y] at points xi.
    Requires x to be sorted in ascending order.
    Extrapolates out of range.
"""
function LinInterp(n, x, y, ni, xi)

    yi      = zeros(ni, 1)
    LocL    = zeros(ni)
    xL = zeros(ni)
    xH = zeros(ni)
    yL = zeros(ni)
    yH = zeros(ni)

    i = 1
    while i<=ni
        if minimum(x) >= xi[i]
            v, LocL[i] = findmax(x)
        else
            v, LocL[i] = findmax(x[xi[i] .> x])
        end
        if xi[i] <= x[1]
            LocL[i] = 1
        end
        if LocL[i] >= n
            LocL[i] = n-1
        end
        xL[i] = x[round(Int, LocL[i])]
        xH[i] = x[round(Int, LocL[i])+1]
        yL[i] = y[round(Int, LocL[i])]
        yH[i] = y[round(Int, LocL[i])+1]
        yi[i] = yL[i] + (xi[i]-xL[i])*((yH[i]-yL[i])/(xH[i]-xL[i]))
        i = i+1
    end

    yi

end

export LinInterp