export FnTaxParamNetMA0NuUI2

using Distributions
using Random

"""
    FnTaxParamNetMA0NuUI2(GrossIncome, Twork, ngpz, ngpe, ngpal, kappa,
        aldist, algrid, edist, egrid, zdist, zgrid,
        stax, ptax, btax, popsize, theta,
        targetTaxToLabinc, Display, UI,
        g0, g1, g2, a, b, c, d)

Build the grid of annual earnings implied by the GKOS (2021) earnings process with a binary nonemployment shock and UI benefits, before and after Gouveia--Strauss taxes, and compute average earnings by age, total labor income and the tax-calibration residual. All arguments are positional.

# Inputs

## Switches and dimensions

- `GrossIncome::Int`: if ``1``, the returned `ygrid` holds pre-tax income; otherwise it holds after-tax income.
- `Twork::Int`: number of working-life periods ``T_{\\text{work}}``.
- `ngpz::Int`: number of grid points for the persistent component ``z``.
- `ngpe::Int`: number of grid points for the transitory component ``\\varepsilon``.
- `ngpal::Int`: number of grid points for the fixed effect ``\\alpha``.
- `Display::Int`: if ``1``, print tax revenue as a percentage of pre-tax labor income.

## Component grids (from the component modules)

- `algrid::Vector` (length ``ngpal``): fixed-effect grid ``\\alpha_i``.
- `aldist::Vector` (length ``ngpal``): probability weights ``\\pi^\\alpha(i)``.
- `zgrid::Matrix` (``Twork`` ``\\times`` ``ngpz``): persistent-component grid ``z_``t,j`` by age.
- `zdist::Matrix` (``Twork`` ``\\times`` ``ngpz``): weights ``\\pi^z_t(j)`` by age.
- `egrid::Matrix` (``Twork`` ``\\times`` ``ngpe``): transitory-component grid ``\\varepsilon_{t,k}`` by age.
- `edist::Matrix` (``Twork`` ``\\times`` ``ngpe``): weights ``\\pi^\\varepsilon_{t(k)}`` by age.
- `popsize::Vector` (length at least ``Twork``): population weight of each age.


## Earnings-process parameters (taken from (Guvenen et al., 2021) and readjusted to account for age rescaling)
- `g0`, `g1`, `g2`: coefficients of the life-cycle
- `a`, `b`, `c`, `d`: coefficients of the nonemployment
- `UI::Real`: unemployment-insurance benefit ``b^{UI}``, in the same units as earnings.

## Tax parameters

- `btax`, `ptax`, `stax`: Gouveia--Strauss parameters ``\\tau_b``, ``\\tau_p``, ``\\tau_s``.
- `targetTaxToLabinc::Real`: target ratio of tax revenue to pre-tax labor income.

## Unused

- `kappa`, `theta`: accepted but not used in the body of the function.


# Outputs

A `{Vector{Any}}` of eight elements, in this order:

- `FnTaxParamNet::Float64`: tax-calibration residual.
- `ypregrid::Array{Float64}`: pre-tax income grid,
- `ygrid::Array{Float64}`: income grid for the household problem
- `pzgrid::Matrix{Float64}`: nonemployment probabilities ``p_\\nu(t, z_{t,j})``,
- `avearnspre::Vector{Float64}`: average pre-tax earnings by age, equation~\\eqref``eq:avpre``.
- `avearnspost::Vector{Float64}`: average after-tax income by age, equation~\\eqref``eq:avpost``.
- `totlabincpre::Float64`: total pre-tax labor income ``\\text{LI^{\\text{pre}}``.
- `totlabincpost::Float64`: total after-tax labor income ``\\text{LI^{\\text{post}}}``.

"""
function FnTaxParamNetMA0NuUI2(GrossIncome, Twork, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, stax, ptax, btax, popsize, theta, targetTaxToLabinc, Display, UI, g0, g1, g2, a, b, c, d)

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
    ygrid = zeros(Twork, ngpz, ngpe, ngpal, 2)
    ypregrid = zeros(Twork, ngpz, ngpe, ngpal, 2)
    ypostgrid = zeros(Twork, ngpz, ngpe, ngpal, 2)

    # gouveia strauss 
    stax = stax;   # guess: is chosen optimally
    ptax = ptax;
    btax = btax;

    #other taxes
    pentax = 0.0;    # payroll tax   
    rtax = 0.0; # tax on interest income
    otax = 0.0; # tax on old age pensions, check option in Parameters
    pencapfrac = 2.2; # cap on [pre-tax] earnings that contribute to pension index, as fraction of av pre-tax earnigs


    g0=g0
    g1=g1
    g2=g2


    a=a
    b=b
    c=c
    d=d

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
                xi=a+b*it+c*zgrid[it, iz]+d*it*zgrid[it, iz]
                xi=max(xi, -200)
                xi=min(xi, 200)
                p=exp(xi)/(1+exp(xi))
                pzgrid[it, iz]=p
                #display(xi)
                #display(p)
                yg = exp(g + algrid[ia] + egrid[it, ie] + zgrid[it, iz])
                if GrossIncome==1
                    ygrid[it, iz, ie, ia, 1] = yg #nu=0 with proba 1-p
                    ygrid[it, iz, ie, ia, 2] = UI
                else
                    ygrid[it, iz, ie, ia, 1] = max(yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax)) + pentax*yg, UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI)
                    ygrid[it, iz, ie, ia, 2] = UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI
                end
                ypregrid[it, iz, ie, ia, 1] = yg
                ypregrid[it, iz, ie, ia, 2] = UI
                ypostgrid[it, iz, ie, ia, 1] = max(yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax)) + pentax*yg, UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI)
                ypostgrid[it, iz, ie, ia, 2] = UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI
                ltotlabincpre = ltotlabincpre + yg*zdist[it, iz]*edist[it, ie]*aldist[ia]*popsize[it]*(1-p)
                ltotlabincpost = ltotlabincpost + ypostgrid[it, iz, ie, ia, 1]*zdist[it, iz]*edist[it, ie]*aldist[ia]*popsize[it]*(1-p) + ypostgrid[it, iz, ie, ia, 2]*zdist[it, iz]*edist[it, ie]*aldist[ia]*popsize[it]*(p)
                #FnTax = btax*(ypregrid[it,iz,ie,ia] - (ypregrid[it,iz,ie,ia]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it,iz,ie,ia]
                #ltottax = ltottax + (FnTax - pentax*ypregrid[it,iz,ie,ia])*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]
                avearnspre[it] = avearnspre[it] + ypregrid[it, iz, ie, ia, 1]*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                avearnspre2[it] = avearnspre2[it] + (ypregrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                avearnspost[it] = avearnspost[it] + ypostgrid[it, iz, ie, ia, 1]*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + ypostgrid[it, iz, ie, ia, 2]*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                avearnspost2[it] = avearnspost2[it] + (ypostgrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + (ypostgrid[it, iz, ie, ia, 2]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                avlearnspre[it] = avlearnspre[it] + log(ypregrid[it, iz, ie, ia, 1])*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                avlearnspost[it] = avlearnspost[it] + log(ypostgrid[it, iz, ie, ia, 1])*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + log(ypostgrid[it, iz, ie, ia, 2])*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                avlearnspost2[it] = avlearnspost2[it] + log(ypostgrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + log(ypostgrid[it, iz, ie, ia, 2]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
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
                    xi=a+b*it+c*zgrid[it, iz]+d*it*zgrid[it, iz]
                    xi=max(xi, -200)
                    xi=min(xi, 200)
                    p=exp(xi)/(1+exp(xi))
                    pzgrid[it, iz]=p
                    #display(xi)
                    #display(p)
                    yg = exp(g + algrid[ia] + egrid[it, ie] + zgrid[it, iz])
                    if GrossIncome==1
                        ygrid[it, iz, ie, ia, 1] = yg #nu=0 with proba 1-p
                        ygrid[it, iz, ie, ia, 2] = UI
                    else

                        ygrid[it, iz, ie, ia, 1] = max(yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax)) + pentax*yg, UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI)
                        ygrid[it, iz, ie, ia, 2] = UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI
                    end
                    ypregrid[it, iz, ie, ia, 1] = yg
                    ypregrid[it, iz, ie, ia, 2] = UI
                    ypostgrid[it, iz, ie, ia, 1] = max(yg - btax*(yg - (yg^(-ptax) + stax)^(-1.0/ptax)) + pentax*yg, UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI)
                    ypostgrid[it, iz, ie, ia, 2] = UI - btax*(UI - (UI^(-ptax) + stax)^(-1.0/ptax)) + pentax*UI
                    ltotlabincpre = ltotlabincpre + yg*zdist[it, iz]*edist[it, ie]*aldist[ia]*popsize[it]*(1-p)
                    ltotlabincpost = ltotlabincpost + ypostgrid[it, iz, ie, ia, 1]*zdist[it, iz]*edist[it, ie]*aldist[ia]*popsize[it]*(1-p) + ypostgrid[it, iz, ie, ia, 2]*zdist[it, iz]*edist[it, ie]*aldist[ia]*popsize[it]*(p)
                    #FnTax = btax*(ypregrid[it,iz,ie,ia] - (ypregrid[it,iz,ie,ia]^(-ptax) + stax)^(-1.0/ptax)) + pentax*ypregrid[it,iz,ie,ia]
                    #ltottax = ltottax + (FnTax - pentax*ypregrid[it,iz,ie,ia])*zdist[it,iz]*edist[it,ie]*aldist[ia]*popsize[it]
                    avearnspre[it] = avearnspre[it] + ypregrid[it, iz, ie, ia, 1]*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                    avearnspre2[it] = avearnspre2[it] + (ypregrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                    avearnspost[it] = avearnspost[it] + ypostgrid[it, iz, ie, ia, 1]*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + ypostgrid[it, iz, ie, ia, 2]*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                    avearnspost2[it] = avearnspost2[it] + (ypostgrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + (ypostgrid[it, iz, ie, ia, 2]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                    avlearnspre[it] = avlearnspre[it] + log(ypregrid[it, iz, ie, ia, 1])*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                    avlearnspre2[it] = avlearnspre2[it] + log(ypregrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p)
                    avlearnspost[it] = avlearnspost[it] + log(ypostgrid[it, iz, ie, ia, 1])*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + log(ypostgrid[it, iz, ie, ia, 2])*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                    avlearnspost2[it] = avlearnspost2[it] + log(ypostgrid[it, iz, ie, ia, 1]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(1-p) + log(ypostgrid[it, iz, ie, ia, 2]^2)*zdist[it, iz]*edist[it, ie]*aldist[ia]*(p)
                    ie=ie+1
                end
                iz=iz+1
            end
            ia=ia+1
        end
        it=it+1
    end

    varearnspre = avearnspre2 - avearnspre .^ 2
    varearnspost = avearnspost2 - avearnspost .^ 2
    varlearnspre = avlearnspre2 - avlearnspre .^ 2
    varlearnspost = avlearnspost2 - avlearnspost .^ 2

    FnTaxParamNet = ltottax/ltotlabincpre - targetTaxToLabinc
    totlabincpre = ltotlabincpre
    totlabincpost = ltotlabincpost

    if Display==1
        display(" Tax revenue / Pre-tax labor income: ")
        display((ltottax/ltotlabincpre)*100)
    end

    [FnTaxParamNet, ypregrid, ygrid, pzgrid, avearnspre, avearnspost, totlabincpre, totlabincpost];
end
