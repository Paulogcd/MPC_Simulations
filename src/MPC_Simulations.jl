module MPC_Simulations

    # Utils
    include("0_Module_FindLinProb1.jl")
    include("0_Module_FnGridAlpha.jl")
    include("0_Module_FnGridKEntrep.jl")
    include("0_Module_FnGridPerm.jl")
    include("0_Module_FnGridTrans.jl")
    include("0_Module_FnSSParam_lowr.jl")
    include("0_Module_FnSSParam.jl")
    include("0_Module_FnTaxParam.jl")
    include("0_Module_golden.jl")
    include("0_Module_LinInterp.jl")
    include("0_Module_LinInterp1.jl")
    include("0_Module_mnbrak.jl")
    include("0_Module_rtnewtGrossInc.jl")
    
    # Require FnSSParam_lowr
    include("0_Module_rtsecSSParam_lowr.jl")
    include("0_Module_rtsecSSParam.jl")

    # Require rtnewtGrossInc
    include("0_Module_FnTaxParamNet.jl")
    include("0_Module_zbrentTaxParamNet.jl")

    # Other
    include("./0_Module_FnTaxParamNetMA0Nu.jl")
    include("./0_Module_FnTaxParamNetMA0NuUI.jl")
    include("./0_Module_FnTaxParamNetMA0NuUI2.jl")

    # Functionized scripts
    include("./1_Function_CashInHands.jl")

    @info("Utils loaded.")
    
end
