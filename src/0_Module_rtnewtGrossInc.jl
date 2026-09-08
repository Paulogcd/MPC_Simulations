"""
Missing documentation
"""
function rtnewtGrossInc(lxguess,lnet,lx1,lx2,xacc,btax,ptax,stax,pentax)

maxit=30
x  = lxguess
x1 = lx1
x2 = lx2

j=1
while j<=maxit
    f = (1.0-pentax - btax)*x + btax*(x^(-ptax) + stax)^(-1.0/ptax) - lnet
    df = (1.0-pentax - btax) + btax* (x^(-ptax-1.0)) * ((x^(-ptax) + stax)^(-1.0/ptax - 1.0))
    dx=f/df
    x=x-dx
    lxguess = x
    if (x1-x)*(x-x2) < 0
        display("rtnewtGrossInc: values jumped out of brackets")
    end
    if abs(dx) < xacc
        return lxguess
    end
    j=j+1
end

display("rtnewt: exceed maximum iterations")
end

export rtnewtGrossInc