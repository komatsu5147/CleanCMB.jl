using CleanCMB
using Healpix
using Optim, LinearAlgebra
using Printf, CSV
using Mmap
using Statistics, Plots
using Tables
# %% Simulation parameters
nrz = 10
ℓmin, ℓmax = 30, 260 # ℓ range for fitting
Alens = 1
Δℓ = 10 # multipole binning size
# %% Foreground cleaning parameters
showβ = true # Show the fitted foreground parameters (for the smoothed covariance)?
β0 = βs0, βd0, Td0 = [-3.0, 1.6, 19.6] # starting foreground parameters for minimisation of -log(likelihood) by `res = optimize(func, [βs0, βd0, Td0])`
ℓswitch = 30 # the multipole below which the foreground parameters are fitted for each band-power
smooth_FWHM = 0.5 # smoothing for the covariance matrix in units of degrees
σβ = σβs, σβd, σTd = [0.5, 0.5, 5] # Gaussian priors: -2*log(likelihood) = (x[1] - βs0)^2 / σβs^2 + (x[2] - βd0)^2 / σβd^2 + (x[3] - Td0)^2 / σTd^2
Σβ = Diagonal(σβ .^ 2)
lnlike_fgprior(x) = -0.5 * (x .- β0)' * (Σβ \ (x .- β0))
# %% Specification of the experiments
ν = [27, 39, 93, 145, 225, 280, 350, 410, 850] # in GHz
νused = [true, true, true, true, true, true, true, true, true] # Which frequencies to use for fitting?
nν = length(ν)
kν = findall(x -> x == true, νused)
FWHM = [91, 63, 30, 17, 11, 9, 0.58, 0.5, 0.23] # in arcmin
σ = FWHM * π / 10800 / √(8 * log(2)) # in radians
uKarcmin = [35, 21, 2.6, 3.3, 6.3, 16, 105, 372, 5.7e5] # in μK arcmin (for temperature; x√2 for pol)
lknee = [30, 30, 50, 50, 70, 100, 700, 700, 700]
αknee = [-2.4, -2.4, -2.5, -3, -3, -3, -1.4, -1.4, -1.4]
# polarisation beam, Eq.(5.8) of Ng & Liu, Int.J.Mod.Phys.D, 8, 61 (1999)
bPl(ℓ, σb) = ifelse(ℓ >= 2, exp(-(ℓ * (ℓ + 1) - 4) * σb^2 / 2), 0)
# %% Read in the hits map and calculate fsky need for the likelihood
nhitsfile = "data/nhits_SAT_r7.FITS"
nhits = readMapFromFITS(nhitsfile, 1, Float64)
nside = nhits.resolution.nside
nhits /= maximum(nhits)
fsky = mean(nhits)^2 / mean(nhits .^ 2)
# %% Read in binned theory power spectra
clt_th = CSV.read("tensor_eebb_binned.csv")
cls_th = CSV.read("scalar_eebb_binned.csv")
ell_eff = clt_th.leff
nbands = length(ell_eff)
rclass = 0.01
# %% Read in covariance matrices and perform parametric maxiimum likelihood fitting
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
    nij = zeros(nν, nν, nbands) # Noise covariance matrix for the ML analysis
    for iν = 1:nν
        ell_fact = ell_eff .* (ell_eff .+ 1) / 2π
        nij[iν, iν, :] =
            2 *
            (uKarcmin[iν] * π / 10800)^2 *
            (1 .+ (ell_eff / lknee[iν]) .^ αknee[iν]) .* ell_fact .*
            bPl.(ell_eff, σ[iν]) .^ -2
    end
    ## Use Parametric Maximum Likelihood method to obtain power spectra of clean maps of the CMB
    A(x) = [ones(length(kν)) synch.(ν[kν], βs = x[1]) dust1.(
        ν[kν],
        βd = x[2],
        Td = x[3],
    )] # Frequency response matrix
    # For ℓ > ℓswitch: Apply parametric maximum likelihood method using a smoothed covariance matrix.
    func_sum(x, cij) =
        -lnlike_fgprior(x) -
        fsky *
        Δℓ *
        sum(
            (2 * ell_eff[jb] + 1) *
            loglike_beta(nij[kν, kν, jb], A(x), cij[kν, kν, jb])
            for jb = 1:nbands
        ) # -log(likelihood) to minimise by `optimize`
    # Smooth covariance matrices to `smooth_FWHM` resolution
    cij = zeros(nν, nν, nbands)
    σsmo = smooth_FWHM * π / 180 / sqrt(8 * log(2))
    for ib = 1:nbands
        cij[:, :, ib] = cov1[:, :, ib] * bPl(ell_eff[ib], σsmo)^2
    end
    res = optimize(x -> func_sum(x, cij), β0)
    if showβ
        println("Fitted parameters: B-mode")
        @show res
    end
    β = Optim.minimizer(res)
    B = [synch.(ν[kν], βs = β[1]) dust1.(ν[kν], βd = β[2], Td = β[3])] # Best-fitting frequency response of the FG
    w_smooth = milca_weights(cov1[kν, kν, :], ones(length(kν)), B)
    cl1[:, irz] = ilc_clean_cij(cov1[kν, kν, :], w_smooth)
    cl2[:, irz] = ilc_clean_cij(cov2[kν, kν, :], w_smooth)
    cl3[:, irz] = ilc_clean_cij(cov3[kν, kν, :], w_smooth)
    # For for ℓ ≤ ℓswitch: Apply parametric maximum likelihood method for each band-power.
    iib = findall(x -> x ≤ ℓswitch, ell_eff)
    if iib ≠ []
        for ib = 1:maximum(iib)
            func(x, cij) =
                -lnlike_fgprior(x) -
                fsky *
                (2 * ell_eff[ib] + 1) *
                loglike_beta(nij[kν, kν, ib], A(x), cij[kν, kν, ib])
            res = optimize(x -> func(x, cov1), β0)
            #@show res
            β = Optim.minimizer(res)
            B = [synch.(ν[kν], βs = β[1]) dust1.(ν[kν], βd = β[2], Td = β[3])]
            w = milca_weights(cov1[kν, kν, ib], ones(length(kν)), B)
            cl1[ib, irz] = ilc_clean_cij(cov1[kν, kν, ib], w)
            cl2[ib, irz] = ilc_clean_cij(cov2[kν, kν, ib], w)
            cl3[ib, irz] = ilc_clean_cij(cov3[kν, kν, ib], w)
        end
    end
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
println("Used νs: ", ν[kν], " GHz")
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
CSV.write("milca_results.csv", t, header = ["irz", "r_wo_FGmarg", "r_w_FGmarg"])
