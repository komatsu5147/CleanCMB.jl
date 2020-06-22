using CSV, Plots
using Statistics
# %% Compare ILC results
ilc = CSV.read("ilc_results_sosatccatp_300sims_seed5147.csv")
println("ILC:")
println(
    "- SAT only: r = ",
    mean(ilc.r_SOonly_wo_FGmarg),
    " ± ",
    std(ilc.r_SOonly_wo_FGmarg),
)
println("- SAT+CCATp: r = ", mean(ilc.r_wo_FGmarg), " ± ", std(ilc.r_wo_FGmarg))
p = scatter(
    ilc.r_wo_FGmarg,
    ilc.r_SOonly_wo_FGmarg,
    xlab = "Tensor-to-scalar ratio, SO-SAT + CCAT-prime",
    ylab = "Tensor-to-scalar ratio, SO-SAT only",
    lab = "",
    legend = :topright,
    xlims = [-0.003, 0.01],
    ylims = [-0.003, 0.01],
    aspect_ratio = :equal,
    title = "Internal Linear Combination",
)
p = plot!(
    ilc.r_wo_FGmarg,
    ilc.r_wo_FGmarg,
    lw = 2,
    lab = "r(SAT only) = r(SAT+CCATp)",
)
p = plot!(
    ilc.r_wo_FGmarg,
    ilc.r_wo_FGmarg .+ 0.93e-3,
    lw = 2,
    lab = "r(SAT only) = r(SAT+CCATp) + 0.93e-3",
)
p = histogram!(
    ilc.r_SOonly_wo_FGmarg,
    inset = (1, bbox(0.5, 0.5, 0.5, 0.4)),
    yticks = nothing,
    subplot = 2,
    bg_inside = nothing,
    lab = "SAT only",
    α = 0.7,
    legend = :best,
    legendfontsize = 5,
    xlims = [-0.005, 0.01],
)
p = histogram!(p[2], ilc.r_wo_FGmarg, α = 0.7, lab = "SAT+CCATp")
p = vline!([0], c = :grey, lw = 2, ls = :dot, lab = "")
p = hline!([0], c = :grey, lw = 2, ls = :dot, lab = "")
savefig("scatterplot_so_vs_soccatp_ilc.pdf")
display(p)

# %% Compare Parametric Maximum Likelihood results
milca = CSV.read("milca_results_sosatccatp_300sims_seed5147.csv")
println("MILCA:")
println(
    "- SAT only: r = ",
    mean(milca.r_SOonly_wo_FGmarg),
    " ± ",
    std(milca.r_SOonly_wo_FGmarg),
)
println(
    "- SAT+CCATp: r = ",
    mean(milca.r_wo_FGmarg),
    " ± ",
    std(milca.r_wo_FGmarg),
)
p = scatter(
    milca.r_wo_FGmarg,
    milca.r_SOonly_wo_FGmarg,
    xlab = "Tensor-to-scalar ratio, SO-SAT + CCAT-prime",
    ylab = "Tensor-to-scalar ratio, SO-SAT only",
    lab = "",
    legend = :topright,
    xlims = [-0.003, 0.01],
    ylims = [-0.003, 0.01],
    aspect_ratio = :equal,
    title = "Parametric Maximum Likelihood",
)
p = plot!(
    milca.r_wo_FGmarg,
    milca.r_wo_FGmarg,
    lw = 2,
    lab = "r(SAT only) = r(SAT+CCATp)",
)
p = plot!(
    milca.r_wo_FGmarg,
    milca.r_wo_FGmarg .+ 0.84e-3,
    lw = 2,
    lab = "r(SAT only) = r(SAT+CCATp) + 0.84e-3",
)
p = histogram!(
    milca.r_SOonly_wo_FGmarg,
    inset = (1, bbox(0.5, 0.5, 0.5, 0.4)),
    yticks = nothing,
    subplot = 2,
    bg_inside = nothing,
    lab = "SAT only",
    α = 0.7,
    legend = :best,
    legendfontsize = 5,
    xlims = [-0.005, 0.01],
)
p = histogram!(p[2], milca.r_wo_FGmarg, α = 0.7, lab = "SAT+CCATp")
p = vline!([0], c = :grey, lw = 2, ls = :dot, lab = "")
p = hline!([0], c = :grey, lw = 2, ls = :dot, lab = "")
savefig("scatterplot_so_vs_soccatp_milca.pdf")
display(p)

# %% Compare ILC and Parametric Maximum Likelihood results
p = scatter(
    ilc.r_wo_FGmarg,
    milca.r_wo_FGmarg,
    xlab = "Tensor-to-scalar ratio, ILC",
    ylab = "Tensor-to-scalar ratio, ML",
    lab = "",
    legend = :topright,
    xlims = [-0.003, 0.01],
    ylims = [-0.003, 0.01],
    aspect_ratio = :equal,
    title = "ILC vs Parametric ML",
)
p = plot!(ilc.r_wo_FGmarg, ilc.r_wo_FGmarg, lw = 2, lab = "r(ML) = r(ILC)")
p = histogram!(
    milca.r_wo_FGmarg,
    inset = (1, bbox(0.5, 0.5, 0.5, 0.4)),
    yticks = nothing,
    subplot = 2,
    bg_inside = nothing,
    lab = "ML",
    α = 0.7,
    legend = :best,
    legendfontsize = 5,
    xlims = [-0.005, 0.01],
)
p = histogram!(p[2], ilc.r_wo_FGmarg, α = 0.7, lab = "ILC")
p = vline!([0], c = :grey, lw = 2, ls = :dot, lab = "")
p = hline!([0], c = :grey, lw = 2, ls = :dot, lab = "")
savefig("scatterplot_ilc_vs_milca.pdf")
display(p)
