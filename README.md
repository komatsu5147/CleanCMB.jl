# CleanCMB

This package contains functions to enable extraction of clean maps of the cosmic microwave background (CMB). Some functions can be used to extract non-CMB astrophysical components as well.

Different algorithms exist for extraction of clean maps of the CMB (as well as of astrophysical components). The package currently supports:

- Internal Linear Combination (ILC) Method
- Parametric Maximum Likelihood Method

See [this note](https://github.com/komatsu5147/CleanCMB.jl/tree/master/note_on_ilc_vs_ml.pdf) for the relationship between them.

## Internal Linear Combination (ILC) Method
- `ilc_weights(cij[, e, ℓid=3])`: return weights for the internal linear combination (ILC) method, following Equation (12) of [Tegmark et al., PRD, 68, 123523 (2003)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.68.123523).
- `cilc_weights(cij, a, b[, ℓid=3])`: return weights for the constrained internal linear combination (CILC) method for two components, following Equation (20) of [Remazeilles, et al., MNRAS, 410, 2481 (2011)](https://academic.oup.com/mnras/article/410/4/2481/1007333).
- `milca_weights(cij, a, B[, ℓid=3])`: return weights for the modified internal linear combination algorithm (MILCA) method (CILC for N components), following Equation (15) of [Hurier, et al., A&A, 558, A118 (2013)](https://www.aanda.org/articles/aa/abs/2013/10/aa21891-13/aa21891-13.html).
  - Note: These papers define weights in various domain, including harmonic, wavelet, and pixel domain. You can choose to work in any domain by supplying a covariance matrix `cij` in the appropriate domain.
- `ilc_clean_cij(cij, w)`: return a covariance matrix multiplied by weights, `w' * cij * w`, with `w` returned by any of the above functions, e.g., `w = cilc_weights(cij, a, b)`.
  - This would be a variance of the extracted component if `cij` were the same as that used in `ilc_weights()`, `cilc_weights()` or `milca_weights()`.
  - If `cij` is a noise covariance matrix, this function returns the noise variance of the cleaned map for each multipole, band-power bin, pixel, etc.

### Arguments
- `cij::Array{<:AbstractFloat,2}`: `nν`-by-`nν` symmetric covariance matrix, where `nν` is the number of frequency bands.
  - or, `cij::Array{<:AbstractFloat,3}`:symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default). Here, `nℓ` is the number of elements in the relevant domain, e.g., multipoles, band-power bins, pixels, etc.
  - The functions `ilc_weights()`, `cilc_weights()` and `milca_weights()` are *multiple dispatch*, which detect the dimension of the input `cij` automatically.
- `a::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. E.g., `a=[1,...,1]` for CMB. The number of elements is `nν`.
- `b::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be nulled. The number of elements is `nν`.
- `B::Array{<:AbstractFloat,2}`: `nν`-by-`nrc` matrix of the frequency response, for `nrc` components to be nulled.

### Optional arguments
- `e::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. The default value is `e=[1,...,1]`, i.e., CMB. The number of elements is `nν`.
- `ℓid::Integer=3`: location of the index for the `nℓ` domain, if the dimension of the input `cij` is `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)`.
`ℓid=1` if `cij[nℓ,nν,nν]`, `ℓid=2` if `cij[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cij[nν,nν,nℓ]`.

## Parametric Maximum Likelihood Method
- `loglike_beta(nij, A, d)` or `loglike_beta(nij, A, cij)`:: return log(likelihood) for frequency response vectors `A` given a data vector `d` or data covariance matrix `cij`, based on Equation (9) of [Stompor et al., MNRAS, 392, 216 (2009)](https://academic.oup.com/mnras/article/392/1/216/1071929).

<!--
- `loglike_beta_deriv(nij, A, dAdβ, d)` and `loglike_beta_hessian(nij, A, dAdβI, dAdβJ, d)`: return the gradient and hessian of log(likelihood) with respect to foreground parameters, using Equation (A1) and (5) of [Errard et al., PRD, 84, 063005 (2011)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.84.063005), respectively. These functions may be used to make the maximisation of `loglike_beta(nij, A, d)` more efficient.
-->

### Arguments
- `nij::Array{<:AbstractFloat,2}`: `nν`-by-`nν` symmetric noise covariance matrix, where `nν` is the number of frequency bands.
- `A::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the frequency response, for `nc` components in sky.
    - E.g., ``A = [a B]`` where `a = ones(nν)` for CMB and `B` is a `nν`-by-`nc-1` matrix for the frequency response of foreground components.
- `d::Array{<:AbstractFloat,1}`: data vector for a given pixel, (ℓ,m), or any other appropriate domain. The number of elements is `nν`.
- `cij::Array{<:AbstractFloat,2}`: `nν`-by-`nν` symmetric covariance matrix for a given multipole, pixel, or any other appropriate domain.

<!--
- `dAdβ::Array{<:AbstractFloat,2}`, `dAdβI::Array{<:AbstractFloat,2}` and `dAdβJ::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the derivative of the frequency response with respect to a foreground parameter.
    - E.g., ``dAdβ = [zeros(nν) dsynch/dβs zeros(nν)]`` where `zeros(nν)` for CMB and dust because they do not depend on the synchrotron index `βs`.
-->

## Foreground models
The package contains the following functions to return frequency dependence of foreground components:
- `tsz(νGHz; units="cmb", Tcmb=2.7255)`: Spectrum of the thermal Sunyaev-Zeldovich effect given in Equation (V) in Appendix of [Zeldovich, Sunyaev, Astrophys. Space Sci. 4, 301 (1969)](http://articles.adsabs.harvard.edu/pdf/1969Ap%26SS...4..301Z).
- `dust1(νGHz; Td=19.6, βd=1.6, νd=353, units="cmb", Tcmb=2.7255)`: Spectrum of 1-component modified black-body thermal dust emission. The output is normalized to unity at `νGHz = νd`.
- `synch(νGHz; βs=-3, νs=23, Cs=0, νC=40, units="cmb", Tcmb=2.7255)`: Spectrum of synchrotron emission. The output is normalized to unity at `νGHz = νs`.

### Arguments
- `νGHz::Real`: frequency in units of GHz.

### Optional keyword arguments
- `units::String`: units of the spectrum. For Rayleigh-Jeans temperature (brightness temperature) units, set `units = "rj"`. The default is the CMB units.
- `Tcmb::Real`: present-day temperature of the CMB in units of Kelvin. The default value is `2.7255`.
- `Td::Real`: dust temperature. The default value is `19.6` (Kelvin).
- `βd::Real`: dust emissivity index. The default value is `1.6`.
- `νd::Real`: frequency at which the output dust spectrum is normalized to unity. The default value is `353` (GHz).
- `βs::Real`: synchrotron power-law index. The default is `-3`.
- `νs::Real`: frequency at which the output synchrotron spectrum is normalized to unity. The default is `23` (GHz).
- `Cs::Real`: curvature of the synchrotron spectrum. The default value is `0`.
- `νC::Real`: pivot frequency for curvature of the synchrotron spectrum. The default value is `40` (GHz).

## Example Julia Codes
- [examples/CleanWMAP.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/CleanWMAP.jl)
  - This code applies `ilc_weights()` in pixel domain to produce a clean map of CMB from the temperature maps of Wilkinson Microwave Anisotropy Probe (WMAP) at five frequency bands.
  - The code requires [Healpix.jl](https://github.com/ziotom78/Healpix.jl).

### Pipelines
Here we provide example codes for *pipelines*, which do everything from generation of simulated maps to estimation of the cosmological parameter (tensor-to-scalar ratio of the primordial gravitational waves) in one go.
- [examples/ILCPipelineSOSAT.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/ILCPipelineSOSAT.jl)
  - This code applies `ilc_weights()` in harmonic domain to produce power spectra of clean polarisation maps of the CMB.
  - The code requires [Healpix.jl](https://github.com/ziotom78/Healpix.jl) and [Libsharp.jl](https://github.com/ziotom78/Libsharp.jl). It also calls python wrappers [pymaster](https://github.com/LSSTDESC/NaMaster) and [classy](https://github.com/lesgourg/class_public/wiki/Python-wrapper) via `PyCall`.
  - This is a simulation pipeline for the Small Aperture Telescope (SAT) of the [Simons Observatory](https://simonsobservatory.org). The code performs the following operations:
    1. Read in a hits map and calculate weights for inhomogeneous noise.
    2. `include("compute_cl_class.jl")` (see [examples/compute_cl_class.jl](https://github.com/komatsu5147/CleanCMB.jl/blob/master/examples/compute_cl_class.jl)) to generate theoretical scalar- and tensor-mode power spectra of the CMB using [CLASS](https://github.com/lesgourg/class_public).
    3. Read in and smooth foreground maps by beams of the telescope.
    4. Generate Gaussian random realisations of the CMB and noise.
    5. Calculate covariance matrices `cij` of the foreground, noise, and CMB+foreground+noise maps.
    6. Calculate ILC weights in harmonic domain using `ilc_weights()`.
    7. Calculate power spectra of the clean CMB maps using `ilc_clean_cij()`.
    8. Show results for visual inspection, if `showresults = true`.
    9. Calculate the tensor-to-scalar ratio and its uncertainty from the simulated realisations.
    10. Write out the results (tensor-to-scalar ratios with and without residual foreground marginalisation) to a file `ilc_results_sosat.csv` and create a PDF figure `ilc_clbb_sim_sosat.pdf` showing the cleaned B-mode power spectrum, minus noisebias, and minus the residual foreground.
  - For your reference, the results from 300 realisations starting with the initial seed of `5147` are given in [examples/results/ilc_results_sosat_300sims_seed5147.csv](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/results/ilc_results_sosat_300sims_seed5147.csv). You can compute the mean and standard deviation of the tensor-to-scalar ratios. You should find, for 30 < ℓ < 260:
    - Without FG marginalisation: r = (2.0 ± 1.6) x 10<sup>-3</sup>
    - With marginalisation: r = (-1.0 ± 2.7) x 10<sup>-3</sup>

- [examples/MILCAPipelineSOSAT.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/MILCAPipelineSOSAT.jl)
  - This code performs a hybrid of the maximum likelihood and ILC methods. See [this note](https://github.com/komatsu5147/CleanCMB.jl/tree/master/note_on_ilc_vs_ml.pdf) for the relationship between them.
    - The code applies `loglike_beta()` in harmonic domain to find the best-fitting synchrotron and dust spectral indices (`βs` and `βd`), and calculates weights using the N-component constrained ILC `milca_weights()`.
    - When calculating the weights, it uses the total covariance matrix `cij` (like for the ILC) rather than the noise covariance matrix `nij` (like for the maximum likelihood). This minimises further the foreground contribution that is not modeled.
  - The code performs the same operations as the above [ILC pipeline](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/ILCPipelineSOSAT.jl) for the Small Aperture Telescope (SAT) of the [Simons Observatory](https://simonsobservatory.org), but replaces the steps f and g with

    6. Calculate the best-fitting `βs` and `βd` by minimising `-loglike_beta()` with respect to them.
    7. Calculate weights using `milca_weights()` with the best-fitting `βs` and `βd`, and calculate power spectra of the clean CMB maps using `ilc_clean_cij()`.

  - More detail of the procedure:
    - The code finds the best-fitting `βs` and `βd` for each band-power at multipoles below the "switching" multipole, `ℓ ≤ ℓswitch` (the default value is `ℓswitch = 50`). This may better handle complex foreground properties on large angular scales.
    - For higher multipoles, the code finds the global `βs` and `βd` using the covariance matrix smoothed to a given resolution specified by `smooth_FWHM` (in units of degrees; the default value is `smooth_FWHM = 3`). This is equivalent to finding `βs` and `βd` on smoothed maps.
  - For your reference, the results from 300 realisations starting with the initial seed of `5147` are given in [examples/results/milca_results_sosat_300sims_seed5147.csv](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/results/milca_results_sosat_300sims_seed5147.csv). You can compute the mean and standard deviation of the tensor-to-scalar ratios. You should find, for 30 < ℓ < 260:
      - Without FG marginalisation: r = (1.7 ± 2.9) x 10<sup>-3</sup>
      - With marginalisation: r = (-0.8 ± 4.3) x 10<sup>-3</sup>

- [examples/ILCPipelineSOSATCCATp.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/ILCPipelineSOSATCCATp.jl)
  - Same as [examples/ILCPipelineSOSAT.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/ILCPipelineSOSAT.jl), but add simulated data of the [CCAT-prime](https://www.ccatobservatory.org) at 350, 410, and 850 GHz with specifications given in Table 1 of [Choi et al., JLTP, 199, 1089 (2020)](https://link.springer.com/article/10.1007/s10909-020-02428-z).

- [examples/MILCAPipelineSOSATCCATp.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examplesMILCAPipelineSOSATCCATp.jl)
  - Same as [examples/MILCAPipelineSOSAT.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/MILCAPipelineSOSAT.jl), but varies the dust temperature `Td` as the third foreground parameter, and use Gaussian priors for `βs`, `βd` and `Td`.
  - Add simulated data of the [CCAT-prime](https://www.ccatobservatory.org) at 350, 410, and 850 GHz with specifications given in Table 1 of [Choi et al., JLTP, 199, 1089 (2020)](https://link.springer.com/article/10.1007/s10909-020-02428-z).

### Splitting Pipelines
The above codes do everything in one go. However, sometimes it is convenient to split the processes. For example, we do not have to re-generate maps and their covariance matrices everytime we make a small modification to the rest of the pipeline.

Here we provide example codes for splitting the pipelines into two pieces: (1) Generation of simulated maps and their covariance matrices, and (2) Application of foreground cleaning methods to the covariance matrices.

- [examples/GenerateCovMatrices.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/GenerateCovMatrices.jl): Perform the steps (a)-(e) of the pipeline and write out the covariance matrices to binary files in arrays of `(nν, nν, nbands)` where `nν` is the number of frequency channels and `nbands` is the number of band-powers. It also writes out the binned scalar and tensor power spectra used to generate the simulations to text files.
- [examples/PerformILC.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/PerformILC.jl) and [examples/PerformMILCA.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/PerformMILCA.jl): Perform the steps (f)-(j) of of the pipeline and write out the estimated tensor-to-scalar ratios to a csv file.
  - **Note**: [examples/PerformMILCA.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/PerformMILCA.jl) varies three parameters, `βs`, `βd` and `Td` with Gaussian priors.

## Acknowledgment

Development of the functions provided in this package was supported in part by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy - EXC-2094 - 390783311 and JSPS KAKENHI Grant Number JP15H05896. The Kavli IPMU is supported by World Premier International Research Center Initiative (WPI), MEXT, Japan.
