using CleanCMB
using Healpix, Libsharp
using PyCall
using Random, Statistics
using Printf
using Plots, LaTeXStrings
using Tables, CSV
# %% Create figures?
showresults = true

# %% Simulation parameters
nrz = 10 # How many realisations?
#Random.seed!(5147) # Initial random number seed. Useful if you need reproducible sequence
rsim = 0 # Tensor-to-scalar ratio used for the simulation
Alens = 1 # Lensing power spectrum amplitude (Alens = 1 for the fiducial)

# %% Specification of the experiment
# Reference: Simons Observatory Collaboration, JCAP, 02, 056 (2019), Table 1.
ν = [27, 39, 93, 145, 225, 280] # in GHz
nν = length(ν)
FWHM = [91, 63, 30, 17, 11, 9] # in arcmin
σ = FWHM * π / 10800 / √(8 * log(2)) # in radians
uKarcmin = [35, 21, 2.6, 3.3, 6.3, 16] # in μK arcmin (for temperature; x√2 for pol)
lknee = [30, 30, 50, 50, 70, 100]
αknee = [-2.4, -2.4, -2.5, -3, -3, -3]

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

# %% Basic parameters for spherical harmonics transform
lmax, mmax = 3 * nside - 1, 3 * nside - 1
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
        bl = exp(-l^2 * σ[iν]^2 / 2)
        for m = 0:l
            ilm = alm_index(alm_info, l, m) + 1
            elm.alm[ilm] *= bl
            blm.alm[ilm] *= bl
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

# %% Compute the scalar and tensor CMB power spectrum
# Call the python wrapper for CLASS, `classy`
# Reference: https://github.com/lesgourg/class_public
classy = pyimport("classy")
# Create an instance of the CLASS wrapper
cosmo = classy.Class()
## Scalar power spectrum
# Create a dictionary of the cosmological parameters
Tcmb = 2.7255e6 # in μK
deg_ncdm = true # 3 massive neutrinos with degenerate mass?
params = Dict(
    "output" => "tCl pCl lCl",
    "modes" => "s",
    "l_max_scalars" => lmax,
    "lensing" => "yes",
    "A_s" => exp(3.0448) * 1e-10,
    "n_s" => 0.96605,
    "k_pivot" => 0.05,
    "h" => 0.6732,
    "omega_b" => 0.022383,
    "omega_cdm" => 0.12011,
    "tau_reio" => 0.0543,
    "N_ncdm" => 1,
)
if deg_ncdm
    push!(params, "m_ncdm" => 0.02, "deg_ncdm" => 3, "N_ur" => 0.00641)
    mν = params["m_ncdm"] * params["deg_ncdm"]
else
    push!(params, "m_ncdm" => 0.06, "N_ur" => 2.0328)
    mν = params["m_ncdm"]
end
# Set the parameters to the cosmological code
cosmo.set(params)
# Run the whole code. Depending on your output, it will call the
# CLASS modules more or less fast. For instance, without any
# output asked, CLASS will only compute background quantities,
# thus running almost instantaneously.
# This is equivalent to the beginning of the `main` routine of CLASS,
# with all the struct_init() methods called.
cosmo.compute()
# Access the lensed scalar cl until l=lmax. It is a dictionnary that contains the fields: ell, tt, te, ee, bb, pp, tp
cls = cosmo.lensed_cl(lmax)
# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
cosmo.empty()
## Tensor power spectrum
rclass = 0.01
params = Dict(
    "output" => "tCl pCl",
    "modes" => "t",
    "l_max_tensors" => lmax,
    "lensing" => "no",
    "r" => rclass,
    "n_t" => 0,
    "A_s" => exp(3.0448) * 1e-10,
    "k_pivot" => 0.05,
    "h" => 0.6732,
    "omega_b" => 0.022383,
    "omega_cdm" => 0.12011,
    "tau_reio" => 0.0543,
    "N_ncdm" => 1,
)
if deg_ncdm
    push!(params, "m_ncdm" => 0.02, "deg_ncdm" => 3, "N_ur" => 0.00641)
    mν = params["m_ncdm"] * params["deg_ncdm"]
else
    push!(params, "m_ncdm" => 0.06, "N_ur" => 2.0328)
    mν = params["m_ncdm"]
end
cosmo.set(params)
cosmo.compute()
clt = cosmo.raw_cl(lmax)

# %% Import NaMaster for the power spectrum analysis on a partial sky
# Reference: Alonso et al.,  MNRAS, 484, 4127 (2019), https://github.com/LSSTDESC/NaMaster
nmt = pyimport("pymaster")
# Read in and apodize mask (C^2 apodization, 10 degrees)
maskfile = "data/mask_apodized_r7.fits"
m = readMapFromFITS(maskfile, 1, Float64)
mask = nmt.mask_apodization(m.pixels, 10.0, apotype = "C2")
# Create binning scheme. We will use 10 multipoles per bandpower.
b = nmt.NmtBin.from_nside_linear(nside, 10, is_Dell = true)
# Array with effective multipole per bandpower
ell_eff = b.get_effective_ells()
nbands = b.get_n_bands()
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
# The function defined below will compute the power spectrum between two
# NmtFields f_a and f_b, using the coupling matrix stored in the
# NmtWorkspace wsp and subtracting the deprojection bias clb.
# Note that the most expensive operations in the MASTER algorithm are
# the computation of the coupling matrix and the deprojection bias. Since
# these two objects are precomputed, this function should be pretty fast!
function compute_master(f_a, f_b, wsp)
    # Compute the power spectrum (a la anafast) of the masked fields
    # Note that we only use n_iter=0 here to speed up the computation,
    # but the default value of 3 is recommended in general.
    cl_coupled = nmt.compute_coupled_cell(f_a, f_b)
    # Decouple power spectrum into bandpowers inverting the coupling matrix
    cl_decoupled = wsp.decouple_cell(cl_coupled)
    return cl_decoupled
end

# %% Loop over realisations
ee1, bb1 = zeros(nbands, nrz), zeros(nbands, nrz) # Cleaned power spectra
ee2, bb2 = zeros(nbands, nrz), zeros(nbands, nrz) # Noise power spectra
ee3, bb3 = zeros(nbands, nrz), zeros(nbands, nrz) # Residual foreground power spectra
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
    f1, f2, f3 = [], [], [] # List of NaMaster fields
    for iν = 1:nν
        for l = 0:lmax
            bl = exp(-l^2 * σ[iν]^2 / 2)
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
                celm.alm[ilm] = elm.alm[ilm] * bl
                cblm.alm[ilm] = blm.alm[ilm] * bl
                nelm.alm[ilm] = √ee * randn(Float64)
                nblm.alm[ilm] = √bb * randn(Float64)
                for m = 1:l
                    ilm = alm_index(alm_info, l, m) + 1
                    celm.alm[ilm] = elm.alm[ilm] * bl
                    cblm.alm[ilm] = blm.alm[ilm] * bl
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
        ell = 0:lmax
        bl = exp.(-ell .^ 2 * σ[iν]^2 / 2)
        ## Create NaMaster fields for the pseudo-Cℓ
        # f1: foreground + CMB + noise
        push!(
            f1,
            nmt.NmtField(
                mask,
                [f_q[iν] + c_q + n_q, f_u[iν] + c_u + n_u],
                purify_b = true,
                beam = bl,
            ),
        )
        # f2: noise
        push!(f2, nmt.NmtField(mask, [n_q, n_u], purify_b = true, beam = bl))
        # f3: foreground
        push!(
            f3,
            nmt.NmtField(mask, [f_q[iν], f_u[iν]], purify_b = true, beam = bl),
        )
    end
    ## Compute covariance matrices of EE and BB
    ce1, cb1 = zeros(nν, nν, nbands), zeros(nν, nν, nbands) # total
    ce2, cb2 = zeros(nν, nν, nbands), zeros(nν, nν, nbands) # noise
    ce3, cb3 = zeros(nν, nν, nbands), zeros(nν, nν, nbands) # foreground
    for iν = 1:nν, jν = iν:nν
        w.compute_coupling_matrix(f1[iν], f1[jν], b)
        cl1 = compute_master(f1[iν], f1[jν], w)
        cl2 = compute_master(f2[iν], f2[jν], w)
        cl3 = compute_master(f3[iν], f3[jν], w)
        ce1[iν, jν, :], cb1[iν, jν, :] = cl1[1, :], cl1[4, :]
        ce2[iν, jν, :], cb2[iν, jν, :] = cl2[1, :], cl2[4, :]
        ce3[iν, jν, :], cb3[iν, jν, :] = cl3[1, :], cl3[4, :]
    end
    ## Perform ILC to obtain power spectra of clean maps of the CMB
    we, wb = ilc_weights(ce1), ilc_weights(cb1)
    ee1[:, irz], bb1[:, irz] = ilc_clean_cij(ce1, we), ilc_clean_cij(cb1, wb)
    ee2[:, irz], bb2[:, irz] = ilc_clean_cij(ce2, we), ilc_clean_cij(cb2, wb)
    ee3[:, irz], bb3[:, irz] = ilc_clean_cij(ce3, we), ilc_clean_cij(cb3, wb)
    # Show power spectra for visual inspection
    if showresults
        p = plot(
            ell_eff,
            cls_th_binned[1, :] .+ (rsim / rclass) * clt_th_binned[1, :],
            title = @sprintf("Realisation # = %03d", irz),
            label = "True EE",
            xaxis = :log,
            yaxis = :log,
            xlab = L"\ell",
            ylab = L"\ell(\ell+1)C_\ell/2\pi",
            ylims = [2e-5, 20],
            xlims = [30, 300],
            legend = :topleft,
        )
        p = plot!(
            ell_eff,
            Alens * cls_th_binned[4, :] .+
            (rsim / rclass) * clt_th_binned[4, :],
            label = "True BB",
        )
        p = plot!(ell_eff, ee1[:, irz], m = :circle, lab = "Cleaned EE")
        p = plot!(ell_eff, bb1[:, irz], m = :circle, lab = "Cleaned BB")
        p = plot!(ell_eff, ee2[:, irz], m = :diamond, lab = "Noise EE")
        p = plot!(ell_eff, bb2[:, irz], m = :diamond, lab = "Noise BB")
        p = plot!(ell_eff, ee3[:, irz], m = :star5, lab = "FG EE")
        p = plot!(ell_eff, bb3[:, irz], m = :star5, lab = "FG BB")
        display(p)
    end
end

# %% Calculate the mean power spectra and variance
ℓmin, ℓmax = 30, 260 # ℓ range for fitting
me1, mb1 = zeros(nbands), zeros(nbands) # Mean total power spectra
ve1, vb1 = zeros(nbands), zeros(nbands) # Variance of total power spectra
me2, mb2 = zeros(nbands), zeros(nbands) # Mean noise power spectra
me3, mb3 = zeros(nbands), zeros(nbands) # Mean residual FG power spectra
for ib = 1:nbands
    me1[ib], mb1[ib] = mean(ee1[ib, :]), mean(bb1[ib, :])
    ve1[ib], vb1[ib] = var(ee1[ib, :]), var(bb1[ib, :])
    me2[ib], mb2[ib] = mean(ee2[ib, :]), mean(bb2[ib, :])
    me3[ib], mb3[ib] = mean(ee3[ib, :]), mean(bb3[ib, :])
end
# Plot and save to "ilc_clbb_sim_sosat.pdf"
if showresults
    ii = findall(x -> x >= ℓmin && x <= ℓmax, ell_eff)
    p = scatter(
        ell_eff[ii],
        mb1[ii],
        m = 5,
        lab = "Cleaned power spectrum",
        xlab = L"\ell",
        ylab = L"\ell(\ell+1)C^{BB}_\ell/2\pi~[\mu K^2]",
        legend = :topleft,
    )
    p = scatter!(ell_eff[ii], mb1[ii] - mb2[ii], m = 5, lab = "- noisebias")
    p = scatter!(
        ell_eff[ii],
        mb1[ii] - mb2[ii] - mb3[ii],
        yerr = sqrt.(vb1[ii]),
        m = 5,
        lab = "- foreground",
    )
    p = plot!(ell_eff[ii], Alens * cls_th_binned[4, ii], lab = "Binned scalar")
    p = plot!(
        ell_eff[ii],
        (rsim / rclass) * clt_th_binned[4, ii],
        lab = "Binned tensor",
        ls = :dash,
    )
    savefig("ilc_clbb_sim_sosat.pdf")
    display(p)
end

# %% Calculate the tensor-to-scalar ratio
r = zeros(nrz)
w = zeros(nrz, 2)
# Joint fit for the tensor-to-scalar ratio and the foreground amplitude
y1 = clt_th_binned[4, ii] / rclass
y2 = mb3[ii]
v = vb1[ii]
Fij = [
    sum(y1 .* y1 ./ v) sum(y1 .* y2 ./ v)
    sum(y2 .* y1 ./ v) sum(y2 .* y2 ./ v)
] # 2x2 Fisher matrix for the tensor-to-scalar ratio and the foreground amplitude
Cij = inv(Fij) # Covariance matrix
for irz = 1:nrz
    x = bb1[ii, irz] .- mb2[ii] .- Alens * cls_th_binned[4, ii]
    r[irz] = sum(x .* y1 ./ v) / sum(y1 .^ 2 ./ v)
    z = [sum(x .* y1 ./ v), sum(x .* y2 ./ v)]
    w[irz, 1:2] = Fij \ z
end
# Report the final results
println("Fitted ℓs: ", ell_eff[ii])
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
CSV.write(
    "ilc_results_sosat.csv",
    t,
    header = ["irz", "r_wo_FGmarg", "r_w_FGmarg"],
)
