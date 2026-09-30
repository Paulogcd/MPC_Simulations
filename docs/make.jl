using Documenter
using MPC_Simulations

makedocs(
    sitename = "MPC_Simulations",
    # format = Documenter.HTML(),
    # modules = [MPC_Simulations]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://github.com/Paulogcd/MPC_Simulations/"
)
