# %% Import NaMaster for the power spectrum analysis on a partial sky
# Reference: Alonso et al.,  MNRAS, 484, 4127 (2019), https://github.com/LSSTDESC/NaMaster
nmt = pyimport("pymaster")
# Read in and apodize mask (C^2 apodization, 10 degrees)
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
