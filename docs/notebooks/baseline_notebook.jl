### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 30aec05b-f715-4860-a87c-e361d219d606
begin 
    using Pkg
    
    Pkg.activate(".")
    
    Pkg.develop(Pkg.PackageSpec(path = "../../."))
    using MPC_Simulations

    Pkg.add("PlutoUI")
    using PlutoUI
    
    for pkg in [
        "Distributions",
        "Roots",
        "StatsFuns",
        "Plots",
        "PlotlyJS",
        "Colors",
        "StatsBase",
        "DataFrames",
        "IterTools",
        "Flatten"]
        Pkg.add(pkg)
    end
    
    using Distributions
    using Roots
    using StatsFuns
    using Plots
    using PlotlyJS
    using Colors
    using StatsBase
    using DataFrames
    using IterTools
    using Flatten
    
    using Distributed
    using DelimitedFiles
    using SharedArrays
    using LinearAlgebra
    using Random
    
    using Dates
    using Base.Threads
    
    Pkg.instantiate()      # installs any missing dependencies
end;

# ╔═╡ 6b16a89b-a71d-45bb-9309-8e34902eb87d
md"""
# Package management
"""

# ╔═╡ 31d477d6-ff4a-4cec-8139-7c852df9600d
PlutoUI.TableOfContents()

# ╔═╡ 223db08c-44e5-417e-8fe2-28b00108ce03
md"""
# Parameters
"""

# ╔═╡ 6d2313f6-1641-475b-a8db-8b00e8d8cd54
@bind a NumberField(default = -2.495)

# ╔═╡ c3a9536a-5842-42fa-bd1f-3cf6f38ea162
@bind b NumberField(default = -0.1037)

# ╔═╡ 46d18f4a-f812-47c4-af5a-20a65d3ab832
@bind c NumberField(default = -5.051)

# ╔═╡ ee998d15-9940-41ed-b7c3-eee3083cd31c
@bind d NumberField(default = -0.1087)

# ╔═╡ 254c72bc-adef-4c8f-ab99-60dbb216e058
@bind RobustTheta NumberField(default = 0)

# ╔═╡ 9826eb12-2c5b-43b2-8635-fedf3200b0bf
@bind SDAlpha NumberField(default = 0.472)

# ╔═╡ 4fde4958-29fe-4898-836d-f17aa7548465
@bind RobustSDTrans1    NumberField(default = 0.762)

# ╔═╡ dc56ca6d-27f0-40c4-8920-ffb82afedc14
@bind RobustSDTrans2    NumberField(default = 0.055)

# ╔═╡ 0b62f49c-51a8-463c-a1b8-3d23e5215b5d
@bind pTrans1   NumberField(default = 0.044)

# ╔═╡ bef37f05-f160-4008-86e5-ddaff5054b3c
@bind ETrans1   NumberField(default = 0.134)

# ╔═╡ a3225123-4828-4748-ab91-dd406a7f8801
@bind ETrans2   NumberField(default = -pTrans1*ETrans1/(1-pTrans1))

# ╔═╡ 366769aa-880e-4edb-8ff9-dc45d02a9058
@bind RobustSDPerm1 NumberField(default = 0.113)

# ╔═╡ 408d6848-a82f-46a1-9cc4-d9d36169862e
@bind RobustSDPerm2 NumberField(default = 0.046)

# ╔═╡ fae24ba3-b204-4070-aae4-2f30f7c38614
@bind pPerm1    NumberField(default = 0.176)

# ╔═╡ f4d13dbb-8bf0-4390-98f7-01c7e6c534db
@bind EPerm1    NumberField(default = -0.524)

# ╔═╡ df89feb6-49a0-4a8c-987f-129d959e2279
@bind EPerm2    NumberField(default = -pPerm1*EPerm1/(1-pPerm1))

# ╔═╡ 4ece9a31-8dde-4087-93fe-1b12ff1eb588
@bind RobustBorrowLimit NumberField(default = - 5.52421885544608/(238.132/218.0555);)

# ╔═╡ e3769858-808a-4593-acad-832c1cc50d1e
@bind UI    NumberField(default = 1.50)

# ╔═╡ a5ee2607-05ff-4396-9ae3-aaf5acb180d3
@bind cbar  NumberField(default = 1.5)

# ╔═╡ 1cf07478-34c5-487e-aa4b-00d35b27a3eb
@bind g0    NumberField(default = 2.74)

# ╔═╡ 717ca8b7-b75a-4ee9-bedd-b6df3922ca05
@bind g1    NumberField(default = 0.062)

# ╔═╡ 6635a224-3b96-4d61-9961-685174e8904d
@bind g2    NumberField(default = -0.0016)

# ╔═╡ b7fb3dd3-fbd9-4783-b834-4a6f3eb81c90
begin
	RunOnLinux = 1;
	Display = 0;
	MatchKYRatio = 0;
	EquilibriumR = 0;
	AnnuityMarkets = 1;
	PreTaxIncome = 0;
	IntitialWealthDist = 0;
	QuadraticPref = 0;
	TaxPensionsGS = 1;
	AgeSpecificVariances = 0;
	Pensions = 1;
	# GrossIncome = 0;
	
	#RANDOM SEED
	
	rng = MersenneTwister(1234);
	
	#ROBUSTNESS CHECKS: to obtain Tables J1 and J2
	RobustDiscountFactor = 1;
	RobustInterestRate = 1.01;
	RobustInitSD = 0.450;
	RobustInitWealth = 1;
	RobustAgeRetir = 29;
	RobustDemographicTrend = 1;
	
	# OPTIONS TO GET BACK TO OLD VERSION
	MatchAggPensionBen = 0;
	UseFinalZPension = 0;
	ScaleBendPoints = 0;
	BendPointsPostTax = 0;
	UseNetIncKY = 0;
	UseNetIncToScaleSS = 0;
	
	# GRIDS DIMENSION - STATE VARIABLES
	ngpal = 4;
	ngpe = 5;
	ngpz = 10;
	ngpa = 50;
    # ngpa = 300;
	ngpm = 5;
	ngpp = ngpal * ngpm * ngpz * ngpe;
	
	# DEMOGRAPHIC PARAMETERS PARAMETERS
	
	Twork = 37;
	Tret = RobustAgeRetir;
	Ttot = Twork+Tret;
	
	# TARGET MOMENTS
	
	targetKY = 2.5;
	targetTaxToLabinc = 0.20;
	targetSSBenToLabinc = 0.095;
	targetSSAvReplacement = 0.45;
	
	# PARAMETERS FOR GRID CONSTRUCTION
	
	pexpgrid = 0.18;
	amax = 300000.0;
	
	# SIMULATION PARAMETERS	
	
	nsim = 50;
	
	# EARNINGS PROCESS 
	
	theta = RobustTheta;
	
	SDz0 = RobustInitSD;
	
	rho = 0.991;
	dPerm = MixtureModel(Normal[
	        Normal(EPerm1, RobustSDPerm1),
	        Normal(EPerm2, RobustSDPerm2)], [pPerm1, 1-pPerm1])
	SDeta = sqrt(var(dPerm));
	
	# INTEREST RATE
	
	R = RobustInterestRate;
	Rsave = RobustInterestRate;
	Rdebt = RobustInterestRate;
	
	# PREFERENCE PARAMETERS 
	
	gam = ones(Ttot);
	
	beta = RobustDiscountFactor;
	betret = RobustDiscountFactor;
	
	qprefb = 120000.0;
	qprefa = 1.0;
	
	# BORROWING LIMIT: SET TO VERY LARGE NEGATIVE NO. FOR NBL
	
	borrowlim = RobustBorrowLimit;
	
	# GOVERNMENT PARAMETERS 
	
	# gouveia strauss 
	# stax = 2.0e-4; # Defined in the next cell
	ptax = 0.768;
	btax = 0.258;
	
	# other taxes and benefits
	cfloor = -100000000000.0;
	pentax = 0.0;
	rtax = 0.0;
	otax = 0.0;
	pencapfrac = 2.2;
	
	Rnet = 1.0 + (1.0-rtax)*(R-1);
	Rnetsave = 1.0 + (1.0-rtax)*(Rsave-1);
	Rnetdebt = 1.0 + (1.0-rtax)*(Rdebt-1);
end;

# ╔═╡ efcd4f9f-f32b-4b68-b973-01622a8155a7
md"""
## Setup initialisation
"""

# ╔═╡ 4d83e431-91ea-4bcd-b7ca-4dcf11a5a992
begin

    # Same-cell definition:
    stax = 2.0e-4
    GrossIncome = 0;
    
    # log experience profile
    
    kappa = open("./input/kappasmoothPSID35.txt") do file
        readlines(file)
    end;
    kappa = map(x->parse(Float64, x), kappa);
    kappa = zeros(Twork);
    
    # demographics profile
    
    bet=ones(Ttot)
    bet[1:Twork]=beta .* ones(Twork)
    bet[(Twork+1):Ttot]=beta .* ones(Tret)
    bet[(Twork-11+1):Ttot]=RobustDemographicTrend*beta .* ones(Tret+11)
    
    # survival probabilities
    
    surprob = open("./input/surprobsmooth.txt") do file
        readlines(file)
    end;
    surprob = map(x->parse(Float64, x), surprob);
    
    # age-specific variances
    
    if AgeSpecificVariances==1
        varetaage = open("./input/varetaage.txt") do file
            readlines(file)
        end;
        varetaage = map(x->parse(Float64, x), varetaage);
    else
        varetaage = ones(Twork-1);
    
    end;
    Vetavec = (SDeta^2)*varetaage;
    
    unconsurprob=zeros(1, Tret);
    annprem=zeros(1, Tret);
    unconsurprob[1] = 1.0;
    
    it=1;
    while it<=Tret
        if it>1
            unconsurprob[it] = unconsurprob[it-1] * surprob[it-1]
        end
        if AnnuityMarkets == 1
            annprem[it] = surprob[it]
        else
            ()
            annprem[it] = 1.0
        end
        global it=it+1
    end;
    
    popsize=zeros(1, Twork+Tret);
    
    it=1;
    while it<=Twork
        popsize[it] = 1.0/real(Twork + sum(unconsurprob))
        global it=it+1
    end;
    
    it=1;
    while it<=Tret
        popsize[Twork+it] = unconsurprob[it]/real(Twork + sum(unconsurprob))
        global it=it+1
    end;
    
    # Initial wealth distribution
    
    initwealthdist = open("./input/initwealthdist.txt") do file
        readdlm(file)
    end;
    
    # Individual fixed effects
    
    lxa = 1.0e-8;
    lxb = 1.0;
    
    lxa, lxb, lxc, lfa, lfb, lfc = mnbrak(lxa, lxb, d -> FnGridAlpha(d, SDAlpha, ngpal)[1]);
    
    lm = golden(lxa, lxb, lxc, 1.0e-2, d -> FnGridAlpha(d, SDAlpha, ngpal)[1]);
    
    lfval, aldist, algrid = FnGridAlpha(lm, SDAlpha, ngpal);
    
    # Transitory shocks
    
    lxa = 1.0e-8;
    lxb = 1.0;
    
    lxa, lxb, lxc, lfa, lfb, lfc = mnbrak(lxa, lxb, d -> FnGridTrans(d, ETrans1, ETrans2, RobustSDTrans1, RobustSDTrans2, pTrans1, ngpe, Twork)[1]);
    
    lm = golden(lxa, lxb, lxc, 1.0e-2, d -> FnGridTrans(d, ETrans1, ETrans2, RobustSDTrans1, RobustSDTrans2, pTrans1, ngpe, Twork)[1]);
    
    lfval, edist, egrid = FnGridTrans(lm, ETrans1, ETrans2, RobustSDTrans1, RobustSDTrans2, pTrans1, ngpe, Twork);
    
    # Permanent component
    
    varz = zeros(Twork);
    varz[1] = SDz0^2;
    
    it=2;
    while it <= Twork
        varz[it] = (rho^2)*varz[it-1] + Vetavec[it-1]
        global it=it+1
    end;
    
    lxa = 1.0e-8;
    lxb = 1.0;
    
    lxa, lxb, lxc, lfa, lfb, lfc = mnbrak(lxa, lxb, d -> FnGridPerm(d, EPerm1, EPerm2, RobustSDPerm1, RobustSDPerm2, pPerm1, Twork, varz, ngpz, rho, SDz0)[1]);
    lm = golden(lxa, lxb, lxc, 1.0e-2, d -> FnGridPerm(d, EPerm1, EPerm2, RobustSDPerm1, RobustSDPerm2, pPerm1, Twork, varz, ngpz, rho, SDz0));
    lfval, varzapprox, zdist, zgrid, ztrans = FnGridPerm(lm, EPerm1, EPerm2, RobustSDPerm1, RobustSDPerm2, pPerm1, Twork, varz, ngpz, rho, SDz0);
    
    # Total income
    
    if PreTaxIncome==0
        GrossIncome, ypregrid, ygrid, pzgrid, avearnspre, avearnspost, totlabincpre, totlabincpost = FnTaxParamNetMA0NuUI2(GrossIncome, Twork, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, stax, ptax, btax, popsize, theta, targetTaxToLabinc, Display, UI, g0, g1, g2, a, b, c, d)
        display(" Found required tax parameter: ")
        display(stax)
    elseif PreTaxIncome==1
        stax, ypregrid, ygrid, avearnspre, avearnspost, totlabincpre, totlabincpost = zbrentTaxParamNet(1.0e-8, 1, 0.000001, d -> FnTaxParam(d, Twork, ngpz, ngpe, kappa, edist, egrid, zdist, zgrid, btax, ptax, pentax, popsize, theta, targetTaxToLabinc, Display))
        display(" Found required tax parameter: ")
        display(stax)
    end;
    
    # Mean pre-tax earnings grid for pension
    
    pencap = pencapfrac*sum(avearnspre)/real(Twork);
    
    # Simulate to get distribution of average pre-tax incomes
    
    zsimI = zeros(Int8, nsim, Twork);
    esimI = zeros(Int8, nsim, Twork);
    alsimI = zeros(Int8, nsim);
    nusimI = zeros(Int8, nsim, Twork);
    
    
    zsimI[:, 1]=sample(rng, collect(1:ngpz), weights(zdist[1, :]), nsim);
    esimI[:, 1]=sample(rng, collect(1:ngpe), weights(edist[1, :]), nsim);
    alsimI[:]=sample(rng, collect(1:ngpal), weights(aldist[:]), nsim);
    
    ypresim = zeros(nsim, Ttot);
    yavsim = zeros(nsim, Twork);
    
    inn=1;
    
    while inn<=nsim
        nudist = Binomial(1, pzgrid[1, zsimI[inn, 1]])
        nu = rand(rng, nudist, 1)[1] #nu=1 means unemployment; nu=0 means employment
        #I compute only on employed people=#
        ypresim[inn, 1] = ypregrid[1, zsimI[inn, 1], esimI[inn, 1], alsimI[inn], 1]
        yavsim[inn, 1] = min(ypresim[inn, 1], pencap)
        it=2
        while it<=Twork
            zsimI[inn, it]=sample(rng, collect(1:ngpz), weights(ztrans[it-1, zsimI[inn, it-1], :]), 1)[1, 1]
            esimI[inn, it]=sample(rng, collect(1:ngpe), weights(edist[it, :]), 1)[1, 1]
            nudist = Binomial(1, pzgrid[it, zsimI[inn, it]])
            nu = rand(rng, nudist, 1)[1]
            #I compute only on employed people=#
            ypresim[inn, it] = ypregrid[it, zsimI[inn, it], esimI[inn, it], alsimI[inn], 1]
            yavsim[inn, it] = ((it-1)*yavsim[inn, it-1] + min(ypresim[inn, it], pencap))/real(it)
            it=it+1
        end
        global inn=inn+1
    end;
    
    # equally spaced between 5th perc and median(), and median and 95th perc
    
    mgrid = zeros(Twork, ngpm);
    
    it=1;
    
    while it<=Twork
        lyavsim = yavsim[:, it]
        lyavsim=sort(lyavsim)
        lperc5 = lyavsim[round(Int, nsim*0.05)]
        lperc50 = lyavsim[round(Int, nsim*0.5)]
        lperc95 = lyavsim[round(Int, nsim*0.95)]
    
        mgrid[it, 1] = lperc5
        mgrid[it, round(Int, (1+ngpm)/2)] = lperc50
        mgrid[it, ngpm] = lperc95
    
        lwidth_1 = (mgrid[it, round(Int, (1+ngpm)/2)] - mgrid[it, 1])/real((ngpm-1)/2)
        lwidth2 = (mgrid[it, ngpm] - mgrid[it, round(Int, (1+ngpm)/2)])/real((ngpm-1)/2)
    
        im=2
        while im<=(ngpm-1)/2
            mgrid[it, im] = mgrid[it, im-1] + lwidth_1
            im=im+1
        end
        im=round(Int, (ngpm+1)/2 + 1)
        while im<=ngpm-1
            mgrid[it, im] = mgrid[it, im-1] + lwidth2
            im=im+1
        end
        global it=it+1
    end;
    
    # Pensions
    
    if Pensions==1
        rtsec, pgrid, ppregrid, pind = rtsecSSParam_lowr(0.1, 1.5, 1.0e-3, Twork, Tret, ngpm, ngpz, ngpe, ngpal, kappa, aldist, algrid, edist, egrid, zdist, zgrid, Display, BendPointsPostTax, ScaleBendPoints, UseFinalZPension, MatchAggPensionBen, avearnspre, avearnspost, mgrid, pencap, targetSSBenToLabinc, targetSSAvReplacement);
    else
        rtsec=1
        pgrid=zeros(Tret, ngpp)
        ppregrid=zeros(Tret, ngpp)
        pind=ones(ngpm, ngpz, ngpe, ngpal)
    end
    
    # Natural borrowing limits
    
    # Last period
    
    nblret = zeros(Tret, ngpp);
    
    ip=1;
    
    while ip<=ngpp
        nblret[Tret, ip] = -pgrid[Tret, ip]/Rnetdebt + 0.0001;
        global ip=ip+1
    end;
    
    # Retired
    
    it=Tret-1;
    
    while it>=1
        ip=1
        while ip<=ngpp
            nblret[it, ip] = annprem[it]*nblret[it+1, ip]/Rnetdebt-pgrid[it, ip]/Rnetdebt + 0.0001;
            ip=ip+1
        end
        global it=it-1
    end;
    
    # Last Working Period
    
    nbl = zeros(Twork, ngpal);
    for ial=1:ngpal
        nbl[Twork, ial] = nblret[1, round(Int, pind[1, 1, 1, ial])]/Rnetdebt - ygrid[Twork, 1, 1, ial, 2]/Rnetdebt + 0.0001;
    end
    
    # Working
    
    it = Twork-1;
    
    while it>=1
        for ial=1:ngpal
            nbl[it, ial] = nbl[it+1, ial]/Rnetdebt - ygrid[it, 1, 1, ial, 2]/Rnetdebt + 0.0001;
        end
        global it=it-1
    end;
    
    # Assets: 0 to amax - exponentially spaced grid() - pexpgrid is parametr
    
    agrid = zeros(Twork, ngpa);
    
    it=1;
    
    while it<=Twork
        agrid[it, 1] = max(borrowlim, nbl[it, 1])
        ia=2
        while ia<=ngpa
            agrid[it, ia] = agrid[it, ia-1] + exp(pexpgrid*(ia-1))
            ia=ia+1
        end
        ltemp1 = (amax-agrid[it, 1])/(agrid[it, ngpa]-agrid[it, 1])
        agrid[it, :] = agrid[it, 1] .+ ltemp1 .* (agrid[it, :] .- agrid[it, 1])
        global it=it+1
    end;
    
    agridret = zeros(Tret, ngpa, ngpp);
    
    it=1;
    
    while it<=Tret
        ip=1
        while ip<=ngpp
            agridret[it, 1, ip] = max(borrowlim, nblret[it, ip])
            ia=2
            while ia<=ngpa
                agridret[it, ia, ip] = agridret[it, ia-1, ip] + exp(pexpgrid*(ia-1))
                ia=ia+1
            end
            ltemp1 = (amax-agridret[it, 1, ip])/(agridret[it, ngpa, ip]-agridret[it, 1, ip])
            agridret[it, :, ip] = agridret[it, 1, ip] .+ ltemp1 .* (agridret[it, :, ip] .- agridret[it, 1, ip])
            ip=ip+1
        end
        global it=it+1
    end;
    
    # Cash in Hand
    
    resultatCashInHands = CashInHands(;Twork, Tret, ngpa, ngpz, ngpe, ngpal, ngpp, ygrid, agrid, agridret, cfloor, Rnetsave, Rnetdebt, pind, pgrid);
    tran        = resultatCashInHands[:tran];
    xgrid       = resultatCashInHands[:xgrid];
    tranret     = resultatCashInHands[:tranret];
    xgridret    = resultatCashInHands[:xgridret];
    
    display("Finished forming grids")
    
    target_mean_a = 238.82286
    facc = 0.02 * target_mean_a
    maxit = 10
    betaL = 0.97343
    betaH = 0.975
    
    j=1
    
    # beta = betaL
    bet = ones(Ttot)
    bet[1:Twork] = betaL .* ones(Twork)
    bet[(Twork-11+1):Ttot] = RobustDemographicTrend * betaL .* ones(Tret+11);
    
end

# ╔═╡ bc62fb6b-b7fd-41c5-a976-20d3c0c33dfd
md"""
## Decisions Simulate30Nucbar_Parallel
"""

# ╔═╡ 8c3afca2-0418-4b95-a9a2-887269d8f8a1
begin 
	resultatAgesResolution = AgesResolution(;Ttot, xgridret, ngpp, qprefa, qprefb, cbar, gam, annprem, pgrid, Display, Tret, QuadraticPref, agridret, Twork, bet, surprob, Rnetdebt, Rnetsave, ngpa, tranret, ngpm, ngpz, ngpe, ngpal, pind, ygrid, xgrid, tran, agrid, mgrid, ypregrid, pencap, pzgrid, ztrans, edist);


	conret      = resultatAgesResolution[:conret];
	mucret      = resultatAgesResolution[:mucret];
	conret1     = resultatAgesResolution[:conret1];
	mucret1     = resultatAgesResolution[:mucret1];
	assret1     = resultatAgesResolution[:assret1];
	muc1        = resultatAgesResolution[:muc1];
	con1        = resultatAgesResolution[:con1];
	ass1        = resultatAgesResolution[:ass1];
	muc         = resultatAgesResolution[:muc];
	con         = resultatAgesResolution[:con];
	ass         = resultatAgesResolution[:ass];
end;

# ╔═╡ 698c5ec7-f2ef-4227-a7f2-cbbfe33c9c89
begin
	Y = con[15, 25, 1, :, 2, 2, 1];
	X = ygrid[15, :, 2, 2, 1];

	Plots.plot(X, Y)
end

# ╔═╡ 4d68a76b-8440-48c2-bb3a-9c4a4d54c84a
md"""
# Simulations
"""

# ╔═╡ 17529f92-010f-42e2-94af-685f2bab86a6
begin
	result = Simulations(nsim, Ttot, rng, ngpal, ngpe, aldist, Twork, edist, ngpz, zdist, ypresim, alsimI, pencap, mgrid, ygrid, initwealthdist, zsimI, esimI, egrid, zgrid, kappa, yavsim, pzgrid, ypregrid, RobustInitWealth, cfloor, Rnetsave, agrid, con, R, it, ztrans, Rnetdebt, pind, agridret, Tret, pgrid, ppregrid, conret, annprem, algrid)

	csim 		= result[:csim]
	csimB 		= result[:csimB]
	csim_mean   = result[:csim_mean]
	csimB_mean 	= result[:csimB_mean]
	csimT_mean  = result[:csimT_mean]
	ysim_mean   = result[:ysim_mean]
	asim_mean   = result[:asim_mean]
	asim 		= result[:asim]
	zsim 		= result[:zsim]
	
end

# ╔═╡ 0d003c38-ca78-4baf-8800-d4dd7e32a1e7
md"""
# Results
"""

# ╔═╡ 29d134e9-1125-4afb-b3b5-e59c3d054ef1
plot_pathCAI(csim_mean, csimT_mean, ysim_mean, asim_mean, Ttot)

# ╔═╡ 0da828e3-c519-42c6-a84d-2276994f1638
md"""
# MPC graphs

1. MPC by wealth (où MPC = la différence entre csimB et csim,  et wealth = asim)

"""

# ╔═╡ 04e3cc38-62c8-4ca3-b280-9bce3ad479ac
begin
	# The mean of the MPC across all individuals
	MPC_toplot_mean = csim_mean-csimB_mean
	# The MPC of all individuals
	MPC_toplot = csim .- csimB
end;

# ╔═╡ 98dff502-edf7-4a82-bc21-27645c11e22e
md"""
A first approach would be to plot the average values, i.e. `MPC_toplot_mean` (y) against `asim_mean` (x). We would obtain:
"""

# ╔═╡ 9d13ec83-0741-478f-81b1-13e464d19155
begin
	Plots.plot(
		asim_mean,
		MPC_toplot_mean,
		xlabel = "wealth",
		ylabel = "MPC",
		label = :none,
		title = "MPC as a function of wealth",
		seriestype = :scatter
	)
end

# ╔═╡ 1be25c52-c7f8-4abb-9a5f-ae24de79a3b0
md"""
We could also plot as a scatter plot the set of individual simulations.

Here, each point corresponds to the situation at one year of one simulated individual.
"""

# ╔═╡ 4934cdc1-2c08-4cd1-bf2e-2574d74ae926
begin
	Plots.plot(
		asim,
		MPC_toplot,
		xlabel = "wealth",
		ylabel = "MPC",
		label = :none,
		title = "MPC as a function of wealth",
		seriestype = :scatter
	)
end

# ╔═╡ 92614c7c-2cdb-48ef-b694-ca03e1a535db
md"""
2. MPC by permanent income (où permanent income = exp(zsim))
"""

# ╔═╡ 1a83505f-22ba-4847-9ae5-c6fc4a52b011
begin
	x_toplot = 1:Ttot
	y_toplot = vcat(mean.(exp.(zsim)[:, i] for i in 1:37), zeros(29))
	z_toplot = MPC_toplot_mean

	matrix_to_plot = hcat(x_toplot, y_toplot, z_toplot)
	df_to_plot = DataFrames.DataFrame(matrix_to_plot, :auto)
	DataFrames.rename!(
		df_to_plot,
		["Age", "Permanent income", "MPC"])
	df_to_plot
end

# ╔═╡ db8d6f55-8804-42f2-b57f-981c3322c6ba
begin
    begin
    n = length(z_toplot)

    plot_MPC_3D = PlotlyJS.Plot(
        PlotlyJS.scatter3d(
            x = x_toplot[1:n],
            y = y_toplot[1:n],
            z = z_toplot[1:n],
            mode = "lines",
            line = attr(width = 4)
        ),
        Layout(
            title = "3D trajectory",
            scene = attr(
                xaxis_title = "Age",
                yaxis_title = "Permanent income",
                zaxis_title = "MPC"
            ),
            width=900, height=600
        )
    )
    plot_MPC_3D
end

end

# ╔═╡ 0ca7c9f1-351f-4200-96c2-dc8b43d5ae03
md"""
# Next steps

- Define clearly the scope of the plots. What do we want to have?
- Functionize the grid creation.
"""

# ╔═╡ Cell order:
# ╟─6b16a89b-a71d-45bb-9309-8e34902eb87d
# ╟─31d477d6-ff4a-4cec-8139-7c852df9600d
# ╟─30aec05b-f715-4860-a87c-e361d219d606
# ╟─223db08c-44e5-417e-8fe2-28b00108ce03
# ╠═6d2313f6-1641-475b-a8db-8b00e8d8cd54
# ╠═c3a9536a-5842-42fa-bd1f-3cf6f38ea162
# ╠═46d18f4a-f812-47c4-af5a-20a65d3ab832
# ╠═ee998d15-9940-41ed-b7c3-eee3083cd31c
# ╠═254c72bc-adef-4c8f-ab99-60dbb216e058
# ╠═9826eb12-2c5b-43b2-8635-fedf3200b0bf
# ╠═4fde4958-29fe-4898-836d-f17aa7548465
# ╠═dc56ca6d-27f0-40c4-8920-ffb82afedc14
# ╠═0b62f49c-51a8-463c-a1b8-3d23e5215b5d
# ╠═bef37f05-f160-4008-86e5-ddaff5054b3c
# ╠═a3225123-4828-4748-ab91-dd406a7f8801
# ╠═366769aa-880e-4edb-8ff9-dc45d02a9058
# ╠═408d6848-a82f-46a1-9cc4-d9d36169862e
# ╠═fae24ba3-b204-4070-aae4-2f30f7c38614
# ╠═f4d13dbb-8bf0-4390-98f7-01c7e6c534db
# ╠═df89feb6-49a0-4a8c-987f-129d959e2279
# ╠═4ece9a31-8dde-4087-93fe-1b12ff1eb588
# ╠═e3769858-808a-4593-acad-832c1cc50d1e
# ╠═a5ee2607-05ff-4396-9ae3-aaf5acb180d3
# ╠═1cf07478-34c5-487e-aa4b-00d35b27a3eb
# ╠═717ca8b7-b75a-4ee9-bedd-b6df3922ca05
# ╠═6635a224-3b96-4d61-9961-685174e8904d
# ╠═b7fb3dd3-fbd9-4783-b834-4a6f3eb81c90
# ╟─efcd4f9f-f32b-4b68-b973-01622a8155a7
# ╠═4d83e431-91ea-4bcd-b7ca-4dcf11a5a992
# ╟─bc62fb6b-b7fd-41c5-a976-20d3c0c33dfd
# ╠═8c3afca2-0418-4b95-a9a2-887269d8f8a1
# ╟─698c5ec7-f2ef-4227-a7f2-cbbfe33c9c89
# ╟─4d68a76b-8440-48c2-bb3a-9c4a4d54c84a
# ╠═17529f92-010f-42e2-94af-685f2bab86a6
# ╟─0d003c38-ca78-4baf-8800-d4dd7e32a1e7
# ╟─29d134e9-1125-4afb-b3b5-e59c3d054ef1
# ╟─0da828e3-c519-42c6-a84d-2276994f1638
# ╠═04e3cc38-62c8-4ca3-b280-9bce3ad479ac
# ╟─98dff502-edf7-4a82-bc21-27645c11e22e
# ╠═9d13ec83-0741-478f-81b1-13e464d19155
# ╟─1be25c52-c7f8-4abb-9a5f-ae24de79a3b0
# ╠═4934cdc1-2c08-4cd1-bf2e-2574d74ae926
# ╟─92614c7c-2cdb-48ef-b694-ca03e1a535db
# ╠═1a83505f-22ba-4847-9ae5-c6fc4a52b011
# ╠═db8d6f55-8804-42f2-b57f-981c3322c6ba
# ╟─0ca7c9f1-351f-4200-96c2-dc8b43d5ae03
