using Plots, Plots.PlotMeasures
using CSV, DataFrames
using Statistics
# %% Read in the fitting results of the Tensor-to-Scalar Ratio (r) of SO-SAT only and SO-SAT+FYST/Prime-Cam
milca_so = CSV.read("sosat/milca_results.csv", DataFrame)
println("MILCA:")
println(
    "- SO-SAT only: r = ",
    mean(milca_so.r_wo_FGmarg),
    " ± ",
    std(milca_so.r_wo_FGmarg),
)
milca = CSV.read("sosatfyst/milca_results.csv", DataFrame)
println(
    "- SO-SAT+FYST: r = ",
    mean(milca.r_wo_FGmarg),
    " ± ",
    std(milca.r_wo_FGmarg),
)
# %% Scatter Plot of r from SO-SAT only and SO-SAT+FYST/Prime-Cam
# Save to "figure7a.pdf"
p = scatter(
    milca.r_wo_FGmarg,
    milca_so.r_wo_FGmarg,
    xlab = "Tensor-to-scalar ratio, SO-SAT + FYST",
    ylab = "Tensor-to-scalar ratio, SO-SAT only",
    lab = "",
    legend = :best,
    xlims = [-0.005, 0.012],
    ylims = [-0.005, 0.012],
    aspect_ratio = :equal,
    markersize = 3,
    legendfontsize = 10,
    tickfontsize = 15,
    labelfontsize = 14,
    #left_margin = -100px,
)
p = plot!(
    milca.r_wo_FGmarg,
    milca.r_wo_FGmarg,
    lw = 2,
    lab = "r(SO-SAT only) = r(SO-SAT+FYST)",
)
p = plot!(
    milca.r_wo_FGmarg,
    milca.r_wo_FGmarg .+ 1.1e-3,
    lw = 2,
    lab = "r(SO-SAT only) = r(SO-SAT+FYST) + 0.0011",
)
savefig("figure7a.pdf")
display(p)
# %% Histogram. Save to "figure7b.pdf"
p = histogram(
    milca_so.r_wo_FGmarg,
    yticks = nothing,
    xlims = [-0.005, 0.012],
    xlab = "Tensor-to-scalar ratio, r",
    lab = "SO-SAT only",
    α = 0.7,
    legend = :best,
    legendfontsize = 15,
    tickfontsize = 13,
    labelfontsize = 20,
)
p = histogram!(milca.r_wo_FGmarg, lab = "SO-SAT+FYST", α = 0.7)
p = vline!([0], c = :grey, lw = 10, lab = "")
savefig("figure7b.pdf")
display(p)
