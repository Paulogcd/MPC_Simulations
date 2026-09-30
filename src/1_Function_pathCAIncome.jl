using Plots
"""
Missing documentation
"""
function plot_pathCAI(csim_mean, csimT_mean, ysim_mean, asim_mean, Ttot)

    result_plot = Plots.plot(
        csim_mean[1:Ttot], label="Path for consumption",
        line=(:blue, 0.5, 6, :solid), size=(800, 600), xticks=(1:(Ttot-1)), ylabel="Consumption",
        xlabel="Year", title="Path for mean consumption",
        xrotation=rad2deg(pi/3), fillrange=0,
        fillalpha=0.25,
        fillcolor=:lightblue, background_color=:ivory
    )

    Plots.plot!(
        result_plot,
        csimT_mean[1:Ttot], label="Path for consumption target",
        line=(:red, 0.5, 6, :solid), size=(800, 600), xticks=(1:(Ttot-1)), ylabel="Consumption",
        xlabel="Year", title="Path for mean consumption",
        xrotation=rad2deg(pi/3), fillrange=0,
        fillalpha=0.25,
        fillcolor=:lightblue, background_color=:ivory
    )

    Plots.plot!(
        result_plot,
        ysim_mean, label="Path for income",
        line=(:green, 0.5, 6, :solid), size=(800, 600), xticks=(1:(Ttot-1)), ylabel="Consumption, Assets, Income",
        xlabel="Year", title="Path for mean consumption, assets, and income",
        xrotation=rad2deg(pi/3), fillrange=0,
        fillalpha=0.25,
        fillcolor=:lightgoldenrod, background_color=:ivory
    )

    Plots.plot!(
        result_plot,
        asim_mean, label="Path for assets",
        line=(:orange, 0.5, 6, :solid), size=(800, 600), xticks=(1:(Ttot-1)), ylabel="Consumption, Assets, Income",
        xlabel="Year", title="Path for mean consumption, assets, and income",
        xrotation=rad2deg(pi/3), fillrange=0,
        fillalpha=0.25,
        fillcolor=:lightgoldenrod, background_color=:ivory
    )
    return result_plot
end
export plot_pathCAI;