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
cosmo.empty()
