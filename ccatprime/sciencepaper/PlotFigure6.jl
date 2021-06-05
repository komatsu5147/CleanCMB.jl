using Plots
using CSV, DataFrames
using Statistics
# %% Read in the fitting results of the Synchrotron Index (βs), Dust Spectral Index (βd) and Dust Temperature (Td)
milca_so = CSV.read("sosat/beta_results.csv", DataFrame)
println("SO-SAT only:")
println("- βs = ", mean(milca_so.beta_s), " ± ", std(milca_so.beta_s))
println("- βd = ", mean(milca_so.beta_d), " ± ", std(milca_so.beta_d))
println("- Td = ", mean(milca_so.T_d), " ± ", std(milca_so.T_d))
milca = CSV.read("sosatfyst/beta_results.csv", DataFrame)
println("SO-SAT+FYST:")
println("- βs = ", mean(milca.beta_s), " ± ", std(milca.beta_s))
println("- βd = ", mean(milca.beta_d), " ± ", std(milca.beta_d))
println("- Td = ", mean(milca.T_d), " ± ", std(milca.T_d))
# %% Scatter Plot of βd and Td. Save to "figure6.pdf"
p = scatter(
    milca_so.beta_d,
    milca_so.T_d,
    xlab = "Dust Spectral Index, βd",
    ylab = "Dust Temperature, Td [K]",
    lab = "SO-SAT only",
    legend = :best,
    xlims = [1.4, 1.9],
    markersize = 3,
    legendfontsize = 15,
    tickfontsize = 15,
    labelfontsize = 18,
)
p = scatter!(milca.beta_d, milca.T_d, lab = "SO-SAT+FYST")
savefig("figure6.pdf")
display(p)
