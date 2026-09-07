"""
Takes in xi, and finds two points in either side of in x and returns the indices of them y and associated probabilities p
"""
function FindLinProb1(x, xi)

    y=zeros(2)
    p=zeros(2)
    n=length(x)
    if minimum(x)>=xi
        v, LocL = findmax(x)
    else
        v, LocL = findmax(x[xi .> x])
    end

    if xi<=x[1]
        y[1]=1
        y[2]=2
        p[1]=1.0
        p[2]=0.0
    elseif LocL>=n
        LocL = n-1
        y[1] = n-1
        y[2] = n
        p[1] = 0.0
        p[2] = 1.0
    else
        y[1] = LocL
        y[2] = LocL+1
        p[1] = (xi - x[LocL])/real(x[LocL+1]-x[LocL])
        p[2] = 1.0 - p[1]
    end

    return [y, p]
end