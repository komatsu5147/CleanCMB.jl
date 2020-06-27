# %% Import NaMaster for the power spectrum analysis on a partial sky
# Reference: Alonso et al.,  MNRAS, 484, 4127 (2019), https://github.com/LSSTDESC/NaMaster
nmt = pyimport("pymaster")
# Read in and apodize mask (C^2 apodization, 10 degrees)
m = readMapFromFITS(maskfile, 1, Float64)
mask = nmt.mask_apodization(m.pixels, 10.0, apotype = "C2")
# Create binning scheme. We will use Δℓ multipoles per bandpower.
b = nmt.NmtBin.from_nside_linear(nside, Δℓ, is_Dell = true)
# Array with effective multipole per bandpower
ell_eff = b.get_effective_ells()
nbands = b.get_n_bands()
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
