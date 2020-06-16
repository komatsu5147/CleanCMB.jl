using CleanCMB
using Printf, CSV
using Mmap
using Statistics, Plots
using Tables
# %% Simulation parameters
nrz = 10
ℓmin, ℓmax = 30, 260 # ℓ range for fitting
Alens = 1
# %% Specification of the experiments
νused = [1, 2, 3, 4, 5, 6, 7, 8, 9] # Which frequencies to use for fitting?
ν = [27, 39, 93, 145, 225, 280, 350, 410, 850] # in GHz
nν = length(ν)
# %% Read in binned theory power spectra
clt_th = CSV.read("tensor_eebb_binned.csv")
cls_th = CSV.read("scalar_eebb_binned.csv")
ell_eff = clt_th.leff
nbands = length(ell_eff)
rclass = 0.01
# %% Read in covariance matrices and perform ILC
io = open("cov_bb_foreg.dat", "r")
cov3 = Mmap.mmap(io, Array{Float64,3}, (nν, nν, nbands))
close(io)
cl1, cl2, cl3 = zeros(nbands, nrz), zeros(nbands, nrz), zeros(nbands, nrz)
for irz = 1:nrz
    io = open(@sprintf("cov_bb_total_irz%03d.dat", irz), "r")
    cov1 = Mmap.mmap(io, Array{Float64,3}, (nν, nν, nbands))
    close(io)
    io = open(@sprintf("cov_bb_noise_irz%03d.dat", irz), "r")
    cov2 = Mmap.mmap(io, Array{Float64,3}, (nν, nν, nbands))
    close(io)
    w = ilc_weights(cov1[νused, νused, :])
    cl1[:, irz] = ilc_clean_cij(cov1[νused, νused, :], w)
    cl2[:, irz] = ilc_clean_cij(cov2[νused, νused, :], w)
    cl3[:, irz] = ilc_clean_cij(cov3[νused, νused, :], w)
end
m1, m2, m3 = mean(cl1, dims = 2), mean(cl2, dims = 2), mean(cl3, dims = 2)
v1 = var(cl1, dims = 2)

# %% Joint fit for the tensor-to-scalar ratio and the foreground amplitude
ii = findall(x -> x >= ℓmin && x <= ℓmax, ell_eff)
r, w = zeros(nrz), zeros(nrz, 2)
y1 = clt_th.bb[ii] / rclass
y2 = m3[ii]
v = v1[ii]
Fij = [
    sum(y1 .* y1 ./ v) sum(y1 .* y2 ./ v)
    sum(y2 .* y1 ./ v) sum(y2 .* y2 ./ v)
] # 2x2 Fisher matrix for the tensor-to-scalar ratio and the foreground amplitude
Cij = inv(Fij) # Covariance matrix
for irz = 1:nrz
    x = cl1[ii, irz] .- m2[ii] .- Alens * cls_th.bb[ii]
    r[irz] = sum(x .* y1 ./ v) / sum(y1 .^ 2 ./ v)
    z = [sum(x .* y1 ./ v), sum(x .* y2 ./ v)]
    w[irz, 1:2] = Fij \ z
end
println("Fitted ℓs: ", ell_eff[ii])
println("Used νs: ", ν[νused], " GHz")
println("Without foreground marginalisation:")
println("- r = ", mean(r), " ± ", std(r))
println("- Fisher error = ", 1 / √sum(y1 .^ 2 ./ v))
println("With foreground marginalisation:")
println("- r = ", mean(w[:, 1]), " ± ", std(w[:, 1]))
println("- Fisher error = ", √Cij[1, 1])
# println("- FG = ", mean(w[:, 2]), " ± ", std(w[:, 2]))
# println("- Fisher error = ", √Cij[2, 2])
# %% Save to the file
t = Tables.table([1:nrz r w[:, 1]])
CSV.write("ilc_results.csv", t, header = ["irz", "r_wo_FGmarg", "r_w_FGmarg"])
