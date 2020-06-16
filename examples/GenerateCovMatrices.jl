using Healpix, Libsharp
using PyCall
using Random, Statistics
using Printf
# %% Simulation parameters
nrz = 10 # How many realisations?
Random.seed!(5147) # Initial random number seed. Useful if you need reproducible sequence
rsim = 0 # Tensor-to-scalar ratio used for the simulation
Alens = 1 # Lensing power spectrum amplitude (Alens = 1 for the fiducial)

# %% Specification of the experiments
# Reference: Simons Observatory Collaboration, JCAP, 02, 056 (2019), Table 1.
# Reference: CCAT-prime Collaboration, J.Low Temp.Phys., 199, 1089 (2020), Table 1.
ν = [27, 39, 93, 145, 225, 280, 350, 410, 850] # in GHz
nν = length(ν)
FWHM = [91, 63, 30, 17, 11, 9, 0.58, 0.5, 0.23] # in arcmin
σ = FWHM * π / 10800 / √(8 * log(2)) # in radians
uKarcmin = [35, 21, 2.6, 3.3, 6.3, 16, 105, 372, 5.7e5] # in μK arcmin (for temperature; x√2 for pol)
lknee = [30, 30, 50, 50, 70, 100, 700, 700, 700]
αknee = [-2.4, -2.4, -2.5, -3, -3, -3, -1.4, -1.4, -1.4]
# temperature beam
bTl(ℓ, σb) = exp(-ℓ * (ℓ + 1) * σb^2 / 2)
# polarisation beam, Eq.(5.8) of Ng & Liu, Int.J.Mod.Phys.D, 8, 61 (1999)
bPl(ℓ, σb) = ifelse(ℓ >= 2, exp(-(ℓ * (ℓ + 1) - 4) * σb^2 / 2), 0)

# %% Read in the hits map and calculate weights for inhomogeneous noise
nhitsfile = "data/nhits_SAT_r7.FITS"
nhits = readMapFromFITS(nhitsfile, 1, Float64)
nside = nhits.resolution.nside
nhits /= maximum(nhits)
fsky = mean(nhits)^2 / mean(nhits .^ 2)
K = √(mean(nhits) / fsky)
weight = zeros(12 * nside^2)
for ip = 1:12*nside^2
    if nhits[ip] > 0
        weight[ip] = K / √nhits[ip]
    end
end

# %% Compute scalar and tensor power spectra using CLASS
lmax, mmax = 3 * nside - 1, 3 * nside - 1
include("compute_cl_class.jl")

# %% Basic parameters for spherical harmonics transform
geom_info = make_healpix_geom_info(nside, 1)
alm_info = make_triangular_alm_info(lmax, mmax, 1)
elm, blm = Alm{ComplexF64}(lmax, mmax), Alm{ComplexF64}(lmax, mmax)

# %% Read in and smooth foreground maps generated by "d1s1" model of Python Sky Model (PySM), in equatorial coordinates
# Reference: Thorne et al., MNRAS, 469, 2821 (2017), https://github.com/healpy/pysm
f_q, f_u = [], []
for iν = 1:nν
    filename = @sprintf("data/map_d1s1_equ_%03dghz_r7_uKcmb_qu.fits", ν[iν])
    q, u = (readMapFromFITS(filename, ic, Float64) for ic = 1:2)
    sharp_execute!(
        SHARP_MAP2ALM,
        2,
        [elm.alm, blm.alm],
        [q.pixels, u.pixels],
        geom_info,
        alm_info,
        SHARP_DP,
    )
    for l = 0:lmax
        for m = 0:l
            ilm = alm_index(alm_info, l, m) + 1
            elm.alm[ilm] *= bPl(l, σ[iν])
            blm.alm[ilm] *= bPl(l, σ[iν])
        end
    end
    sharp_execute!(
        SHARP_ALM2MAP,
        2,
        [elm.alm, blm.alm],
        [q.pixels, u.pixels],
        geom_info,
        alm_info,
        SHARP_DP,
    )
    push!(f_q, q)
    push!(f_u, u)
end

# %% Setup NaMaster for the power spectrum analysis on a partial sky
maskfile = "data/mask_apodized_r7.fits"
include("setup_namaster.jl")
# The theory power spectrum must be binned into bandpowers in the same manner the data has.
# Generate an NmtWorkspace object that we use to compute and store the mode coupling matrix.
# Note that this matrix depends only on the masks of the two fields to correlate, but not on the maps themselves.
w = nmt.NmtWorkspace()
f = nmt.NmtField(mask, [f_q[1], f_u[1]], purify_b = true)
w.compute_coupling_matrix(f, f, b)
cls_th = [cls["ee"], zero(cls["ee"]), zero(cls["bb"]), cls["bb"]] * Tcmb^2
cls_th_binned = w.decouple_cell(w.couple_cell(cls_th))
clt_th = [clt["ee"], zero(clt["ee"]), zero(clt["bb"]), clt["bb"]] * Tcmb^2
clt_th_binned = w.decouple_cell(w.couple_cell(clt_th))
CSV.write(
    "scalar_eebb_binned.csv",
    Tables.table([ell_eff cls_th_binned[1, :] cls_th_binned[4, :]]),
    header = ["leff", "ee", "bb"],
)
clt_th = [clet, zero(clet), zero(clbt), clbt]
clt_th_binned = w.decouple_cell(w.couple_cell(clt_th))
CSV.write(
    "tensor_eebb_binned.csv",
    Tables.table([ell_eff clt_th_binned[1, :] clt_th_binned[4, :]]),
    header = ["leff", "ee", "bb"],
)

# %% Compute covariance matrices of the foreground EE and BB
# Used to compute the power spectra of residual foreground for checking the results,
# but this information cannot be used for the real data analysis.
f3 = []
for iν = 1:nν
    push!(
        f3,
        nmt.NmtField(
            mask,
            [f_q[iν], f_u[iν]],
            purify_b = true,
            beam = bPl.(0:lmax, σ[iν]),
        ),
    )
end
ce3, cb3 = zeros(nν, nν, nbands), zeros(nν, nν, nbands) # foreground
for iν = 1:nν, jν = iν:nν
    w.compute_coupling_matrix(f3[iν], f3[jν], b)
    cl3 = compute_master(f3[iν], f3[jν], w)
    ce3[iν, jν, :], cb3[iν, jν, :] = cl3[1, :], cl3[4, :]
end
# Write out to binary files
open("cov_bb_foreg.dat", "w") do io
    write(io, cb3)
end
open("cov_ee_foreg.dat", "w") do io
    write(io, ce3)
end

# %% Loop over realisations
for irz = 1:nrz
    @show irz
    # Create spherical harmonics coefficients of the CMB
    for l = 2:lmax
        ee = cls_th[1][l+1] + (rsim / rclass) * clt_th[1][l+1]
        bb = Alens * cls_th[4][l+1] + (rsim / rclass) * clt_th[4][l+1]
        ilm = alm_index(alm_info, l, 0) + 1
        elm.alm[ilm] = √ee * randn(Float64)
        blm.alm[ilm] = √bb * randn(Float64)
        for m = 1:l
            ilm = alm_index(alm_info, l, m) + 1
            elm.alm[ilm] = √ee * randn(ComplexF64)
            blm.alm[ilm] = √bb * randn(ComplexF64)
        end
    end
    # Smooth CMB and create spherical harmonics coefficients of the noise
    celm, cblm = Alm{ComplexF64}(lmax, mmax), Alm{ComplexF64}(lmax, mmax)
    nelm, nblm = Alm{ComplexF64}(lmax, mmax), Alm{ComplexF64}(lmax, mmax)
    c_q, c_u = Map{Float64,RingOrder}(nside), Map{Float64,RingOrder}(nside)
    n_q, n_u = Map{Float64,RingOrder}(nside), Map{Float64,RingOrder}(nside)
    f1, f2 = [], [] # List of NaMaster fields
    for iν = 1:nν
        for l = 0:lmax
            ee =
                2 *
                (uKarcmin[iν] * π / 10800)^2 *
                (1 + (l / lknee[iν])^αknee[iν])
            bb =
                2 *
                (uKarcmin[iν] * π / 10800)^2 *
                (1 + (l / lknee[iν])^αknee[iν])
            if l >= 2
                ilm = alm_index(alm_info, l, 0) + 1
                celm.alm[ilm] = elm.alm[ilm] * bPl(l, σ[iν])
                cblm.alm[ilm] = blm.alm[ilm] * bPl(l, σ[iν])
                nelm.alm[ilm] = √ee * randn(Float64)
                nblm.alm[ilm] = √bb * randn(Float64)
                for m = 1:l
                    ilm = alm_index(alm_info, l, m) + 1
                    celm.alm[ilm] = elm.alm[ilm] * bPl(l, σ[iν])
                    cblm.alm[ilm] = blm.alm[ilm] * bPl(l, σ[iν])
                    nelm.alm[ilm] = √ee * randn(ComplexF64)
                    nblm.alm[ilm] = √bb * randn(ComplexF64)
                end
            else
                for m = 0:l
                    ilm = alm_index(alm_info, l, m) + 1
                    celm.alm[ilm], cblm.alm[ilm] = 0, 0
                    nelm.alm[ilm], nblm.alm[ilm] = 0, 0
                end
            end
        end
        sharp_execute!(
            SHARP_ALM2MAP,
            2,
            [celm.alm, cblm.alm],
            [c_q.pixels, c_u.pixels],
            geom_info,
            alm_info,
            SHARP_DP,
        )
        sharp_execute!(
            SHARP_ALM2MAP,
            2,
            [nelm.alm, nblm.alm],
            [n_q.pixels, n_u.pixels],
            geom_info,
            alm_info,
            SHARP_DP,
        )
        n_q .*= weight
        n_u .*= weight
        ## Create NaMaster fields for the pseudo-Cℓ
        # f1: foreground + CMB + noise
        push!(
            f1,
            nmt.NmtField(
                mask,
                [f_q[iν] + c_q + n_q, f_u[iν] + c_u + n_u],
                purify_b = true,
                beam = bPl.(0:lmax, σ[iν]),
            ),
        )
        # f2: noise
        push!(
            f2,
            nmt.NmtField(
                mask,
                [n_q, n_u],
                purify_b = true,
                beam = bPl.(0:lmax, σ[iν]),
            ),
        )
    end
    ## Compute covariance matrices of EE and BB
    ce1, cb1 = zeros(nν, nν, nbands), zeros(nν, nν, nbands) # total
    ce2, cb2 = zeros(nν, nν, nbands), zeros(nν, nν, nbands) # noise
    for iν = 1:nν, jν = iν:nν
        w.compute_coupling_matrix(f1[iν], f1[jν], b)
        cl1 = compute_master(f1[iν], f1[jν], w)
        cl2 = compute_master(f2[iν], f2[jν], w)
        ce1[iν, jν, :], cb1[iν, jν, :] = cl1[1, :], cl1[4, :]
        ce2[iν, jν, :], cb2[iν, jν, :] = cl2[1, :], cl2[4, :]
    end
    # Write out to binary files
    open(@sprintf("cov_bb_total_irz%03d.dat", irz), "w") do io
        write(io, cb1)
    end
    open(@sprintf("cov_bb_noise_irz%03d.dat", irz), "w") do io
        write(io, cb2)
    end
    open(@sprintf("cov_ee_total_irz%03d.dat", irz), "w") do io
        write(io, ce1)
    end
    open(@sprintf("cov_ee_noise_irz%03d.dat", irz), "w") do io
        write(io, ce2)
    end
end
