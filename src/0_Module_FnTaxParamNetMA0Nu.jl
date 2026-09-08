export FnTaxParamNetMA0Nu

using Distributions 
using Random

"""
Documentation missing
"""
function FnTaxParamNetMA0Nu(GrossIncome,Twork,ngpz,ngpe,ngpal,kappa,aldist,algrid,edist,egrid,zdist,zgrid,btax,ptax,pentax,popsize,theta,targetTaxToLabinc,Display)
    
#RANDOM SEED

rng = MersenneTwister(1234);

avearnspre = zeros(Twork)
avearnspre2 = zeros(Twork)
avearnspost = zeros(Twork)
avearnspost2 = zeros(Twork)
avlearnspre = zeros(Twork)
avlearnspre2 = zeros(Twork)
avlearnspost = zeros(Twork)
avlearnspost2 = zeros(Twork)
pzgrid = zeros(Twork, ngpz)
ygrid = zeros(Twork, ngpz, ngpe, ngpal,2)
ypregrid = zeros(Twork, ngpz, ngpe, ngpal,2)
ypostgrid = zeros(Twork, ngpz, ngpe, ngpal,2)
    
# gouveia strauss 
stax = 2.0 #2.0e-4;   					# guess: is chosen optimally
ptax = 0.768;
btax = 0.258;
    
#=quadratic time trend (I divide the estimates that multiply it by 10 because GFOS normalize t=(t-24)/10 while I just have t=t-24)=#
g0=2.746
g1=0.0624
g2=-0.00167

#=proba of nonemployment depends on z and t (I divide the estimates that multiply it by 10 because GFOS normalize t=(t-24)/10 while I just have t=t-24)=#
a=-2.495
b=-0.1037
c=-5.051
d=-0.1087

ltottax = 0.0
ltotlabincpre = 0.0
ltotlabincpost = 0.0

it=1
    avearnspre[it] = 0.0
    avearnspre2[it] = 0.0
    avearnspost[it] = 0.0
    avearnspost2[it] = 0.0
    avlearnspre[it] = 0.0
    avlearnspre2[it] = 0.0
    avlearnspost[it] = 0.0
    avlearnspost2[it] = 0.0
    ia=1
    while ia<=ngpal
    iz=1
    while iz<=ngpz
        ie=1
        while ie<=ngpe
            iep=1
            g=g0+g1*it+g2*it*it
            xi=a+b*it+c*zgrid[it,iz]+d*it*zgrid[it,iz]
            xi=max(xi,-200)
            xi=min(xi,200)
            p=exp(xi)/(1+exp(xi))
            pzgrid[it,iz]=p
            #display(xi)
            display(p)
            yg = exp(g + algrid[ia] + egrid[it,ie] + zgrid[it,iz])
            if GrossIncome==1
                ygrid[it,iz,ie,ia,1] = yg #nu=0 with proba 1-p
                ygrid[it,iz,ie,ia,2] = 0
            else
                
                ygrid[it,iz,ie,ia,1] = yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax))
                ygrid[it,iz,ie,ia,2] = 0 - btax*(0 - (0^(-ptax) + stax)^(-1.0/ptax))
            end
            ypregrid[it,iz,ie,ia,1] = yg
            ypregrid[it,iz,ie,ia,2] = 0
            ypostgrid[it,iz,ie,ia,1] = yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax))
            ypostgrid[it,iz,ie,ia,2] = - btax*(0 - (0^(-ptax) + stax)^(-1.0/ptax))
            ltotlabincpre = ltotlabincpre + yg*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]*(1-p)
            ltotlabincpost = ltotlabincpost + ypostgrid[it,iz,ie,ia,1]*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]*(1-p) + ypostgrid[it,iz,ie,ia,2]*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]*(p)
            FnTax = 0 #btax*(ypregrid[it,iz,ie,ia] - (ypregrid[it,iz,ie,ia]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it,iz,ie,ia]
            ltottax = 0 #ltottax + (FnTax - pentax*ypregrid[it,iz,ie,ia])*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]
            avearnspre[it] = avearnspre[it] + ypregrid[it,iz,ie,ia,1]*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avearnspre2[it] = avearnspre2[it] + (ypregrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avearnspost[it] = avearnspost[it] + ypostgrid[it,iz,ie,ia,1]*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + ypostgrid[it,iz,ie,ia,2]*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            avearnspost2[it] = avearnspost2[it] + (ypostgrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + (ypostgrid[it,iz,ie,ia,2]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            avlearnspre[it] = avlearnspre[it] + log(ypregrid[it,iz,ie,ia,1])*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avlearnspost[it] = avlearnspost[it] + log(ypostgrid[it,iz,ie,ia,1])*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + log(ypostgrid[it,iz,ie,ia,2])*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            avlearnspost2[it] = avlearnspost2[it] + log(ypostgrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + log(ypostgrid[it,iz,ie,ia,2]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            ie=ie+1
        end
        iz=iz+1
    end
    ia=ia+1
    end
it=2

while it<=Twork
    avearnspre[it] = 0.0
    avearnspre2[it] = 0.0
    avearnspost[it] = 0.0
    avearnspost2[it] = 0.0
    avlearnspre[it] = 0.0
    avlearnspre2[it] = 0.0
    avlearnspost[it] = 0.0
    avlearnspost2[it] = 0.0
    ia=1
    while ia<=ngpal
    iz=1
    while iz<=ngpz
        ie=1
        while ie<=ngpe
            g=g0+g1*it+g2*it*it
            xi=a+b*it+c*zgrid[it,iz]+d*it*zgrid[it,iz]
            xi=max(xi,-200)
            xi=min(xi,200)
            p=exp(xi)/(1+exp(xi))
            pzgrid[it,iz]=p
            #display(xi)
            #display(p)
            yg = exp(g + algrid[ia] + egrid[it,ie] + zgrid[it,iz])
            if GrossIncome==1
                ygrid[it,iz,ie,ia,1] = yg #nu=0 with proba 1-p
                ygrid[it,iz,ie,ia,2] = 0
            else
                
                ygrid[it,iz,ie,ia,1] = yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax))
                ygrid[it,iz,ie,ia,2] = 0 - btax*(0 - (0^(-ptax) + stax)^(-1.0/ptax))
            end
            ypregrid[it,iz,ie,ia,1] = yg
            ypregrid[it,iz,ie,ia,2] = 0
            ypostgrid[it,iz,ie,ia,1] = yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax))
            ypostgrid[it,iz,ie,ia,2] = - btax*(0 - (0^(-ptax) + stax)^(-1.0/ptax))
            ltotlabincpre = ltotlabincpre + yg*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]*(1-p)
            ltotlabincpost = ltotlabincpost + ypostgrid[it,iz,ie,ia,1]*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]*(1-p) + ypostgrid[it,iz,ie,ia,2]*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]*(p)
            FnTax = 0 #btax*(ypregrid[it,iz,ie,ia] - (ypregrid[it,iz,ie,ia]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it,iz,ie,ia]
            ltottax = 0 #ltottax + (FnTax - pentax*ypregrid[it,iz,ie,ia])*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]
            avearnspre[it] = avearnspre[it] + ypregrid[it,iz,ie,ia,1]*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avearnspre2[it] = avearnspre2[it] + (ypregrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avearnspost[it] = avearnspost[it] + ypostgrid[it,iz,ie,ia,1]*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + ypostgrid[it,iz,ie,ia,2]*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            avearnspost2[it] = avearnspost2[it] + (ypostgrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + (ypostgrid[it,iz,ie,ia,2]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            avlearnspre[it] = avlearnspre[it] + log(ypregrid[it,iz,ie,ia,1])*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p)
            avlearnspost[it] = avlearnspost[it] + log(ypostgrid[it,iz,ie,ia,1])*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + log(ypostgrid[it,iz,ie,ia,2])*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            avlearnspost2[it] = avlearnspost2[it] + log(ypostgrid[it,iz,ie,ia,1]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(1-p) + log(ypostgrid[it,iz,ie,ia,2]^2)*zdist[it,iz]*edist[it,ie]*aldist[ia]*(p)
            ie=ie+1
        end
        iz=iz+1
    end
    ia=ia+1
    end
    it=it+1
end

varearnspre = avearnspre2 - avearnspre.^2
varearnspost = avearnspost2 - avearnspost.^2
varlearnspre = avlearnspre2 - avlearnspre.^2
varlearnspost = avlearnspost2 - avlearnspost.^2

FnTaxParamNet = ltottax/ltotlabincpre - targetTaxToLabinc
totlabincpre = ltotlabincpre
totlabincpost = ltotlabincpost

if Display==1
    display(" Tax revenue / Pre-tax labor income: ")
    display((ltottax/ltotlabincpre)*100)
    end
        
[FnTaxParamNet,ypregrid, ygrid, pzgrid, avearnspre,avearnspost,totlabincpre,totlabincpost];
end