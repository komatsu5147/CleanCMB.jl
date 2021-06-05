using CleanCMB
using Healpix
using Optim, LinearAlgebra
using Printf, CSV, DataFrames
using Mmap
using Statistics
using Tables
dir = "../data"
# %% Simulation parameters
nrz = 1500
ℓmin, ℓmax = 30, 260 # ℓ range for fitting
ℓko = 100
Alens = 1
Δℓ = 10 # multipole binning size
# %% Foreground cleaning parameters
showβ = false # Show the fitted foreground parameters (for the smoothed covariance)?
β0 = βs0, βd0, Td0 = [-3.0, 1.6, 19.6] # starting foreground parameters for minimisation of -log(likelihood) by `res = optimize(func, [βs0, βd0, Td0])`
ℓswitch = 30 # the multipole below which the foreground parameters are fitted for each band-power
smooth_FWHM = 0.5 # smoothing for the covariance matrix in units of degrees
σβ = σβs, σβd, σTd = [0.5, 0.5, 5] # Gaussian priors: -2*log(likelihood) = (x[1] - βs0)^2 / σβs^2 + (x[2] - βd0)^2 / σβd^2 + (x[3] - Td0)^2 / σTd^2
Σβ = Diagonal(σβ .^ 2)
lnlike_fgprior(x) = -0.5 * (x .- β0)' * (Σβ \ (x .- β0))
# %% Specification of the experiments
ν = [27, 39, 93, 145, 225, 280, 350, 410, 850] # in GHz
νused = [true, true, true, true, true, true, true, true, true] # Which frequencies to use for fitting?
#νused = [true, true, true, true, true, true, false, false, false]
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
nhitsfile = joinpath(dir, "nhits_SAT_r7.FITS")
nhits = readMapFromFITS(nhitsfile, 1, Float64)
nside = nhits.resolution.nside
nhits /= maximum(nhits)
fsky = mean(nhits)^2 / mean(nhits .^ 2)
# %% Read in binned theory power spectra
clt_th = CSV.read(joinpath(dir, "tensor_eebb_binned.csv"), DataFrame)
cls_th = CSV.read(joinpath(dir, "scalar_eebb_binned.csv"), DataFrame)
ell_eff = clt_th.leff
nbands = length(ell_eff)
rclass = 0.01
# %% Read in covariance matrices and perform parametric maxiimum likelihood fitting
io = open(joinpath(dir, "cov_bb_foreg.dat"), "r")
cov3 = Mmap.mmap(io, Array{Float64,3}, (nν, nν, nbands))
close(io)
βout = zeros(3, nrz)
for irz = 1:nrz
    ioin = open(joinpath(dir, @sprintf("cov_bb_total_irz%04d.dat", irz)), "r")
    cov1 = Mmap.mmap(ioin, Array{Float64,3}, (nν, nν, nbands))
    close(ioin)
    ioin = open(joinpath(dir, @sprintf("cov_bb_noise_irz%04d.dat", irz)), "r")
    cov2 = Mmap.mmap(ioin, Array{Float64,3}, (nν, nν, nbands))
    close(ioin)
    nij = zeros(nν, nν, nbands) # Noise covariance matrix for the ML analysis
    for iν = 1:nν
        ell_fact = ell_eff .* (ell_eff .+ 1) / 2π
        nij[iν, iν, :] =
            2 *
            (uKarcmin[iν] * π / 10800)^2 *
            (1 .+ (ell_eff / lknee[iν]) .^ αknee[iν]) .* ell_fact .*
            bPl.(ell_eff, σ[iν]) .^ -2
    end
    ## Knock out low multipole data of FYST
    iko = findall(x -> x <= ℓko, ell_eff)
    for iν = 7:9
        nij[iν, iν, iko] .= Inf
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
            loglike_beta(nij[kν, kν, jb], A(x), cij[kν, kν, jb]) for
            jb = 1:nbands
        ) # -log(likelihood) to minimise by `optimize`
    # Smooth covariance matrices to `smooth_FWHM` resolution
    cij = zeros(nν, nν, nbands)
    cov1ko = copy(cov1)
    cov2ko = copy(cov2)
    σsmo = smooth_FWHM * π / 180 / sqrt(8 * log(2))
    ## Knock out low multipole data of FYST
    for iν = 7:9
        cov1ko[iν, iν, iko] =
            (cov1[iν, iν, iko] .- cov2[iν, iν, iko]) .+ 1e4 * cov2[iν, iν, iko]
        cov2ko[iν, iν, iko] .*= 1e4
    end
    for ib = 1:nbands
        cij[:, :, ib] = cov1ko[:, :, ib] * bPl(ell_eff[ib], σsmo)^2
    end
    res = optimize(x -> func_sum(x, cij), β0)
    if showβ
        println("Fitted parameters: B-mode")
        @show res
    end
    β = Optim.minimizer(res)
    βout[1:3, irz] = β[1:3]
end
# %% Save to the file
t = Tables.table([1:nrz βout[1,:] βout[2,:] βout[3,:]])
CSV.write("beta_results_ko.csv", t, header = ["irz", "beta_s", "beta_d", "T_d"])
