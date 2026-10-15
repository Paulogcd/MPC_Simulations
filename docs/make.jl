using Documenter
using MPC_Simulations


# Here, we could include a function that runs the baseline notebook
# each time the package is updated. We have to make if after the docs
# because makedocs() cleans everything in ./build.

using PlutoSliderServer
# Baseline notebook:
baseline_notebook_path = joinpath(@__DIR__, "notebooks", "baseline_notebook.jl")
PlutoSliderServer.export_notebook(
    baseline_notebook_path;
    Export_output_dir = joinpath(@__DIR__, "src"))

# For debugging purposes:
@show pwd()
@show readdir(joinpath(@__DIR__, "src"))
@show isfile(joinpath(@__DIR__, "src", "baseline_notebook.html"))


Documenter.makedocs(
    sitename = "MPC_Simulations.jl",
    # format = Documenter.HTML(),
    # modules = [MPC_Simulations]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
Documenter.deploydocs(
    repo = "https://github.com/Paulogcd/MPC_Simulations.jl/"
)
