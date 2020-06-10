# CleanCMB

This package contains functions to enable extraction of clean maps of the cosmic microwave background (CMB). Some functions can be used to extract non-CMB astrophysical components as well.

Different algorithms exist for extraction of clean maps of the CMB (as well as of astrophysical components). The package currently supports the following algorithms:

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
- `loglike_beta(nij, A, d)`: return log(likelihood) for frequency response vectors `A` given a data vector `d`, using Equation (9) of [Stompor et al., MNRAS, 392, 216 (2009)](https://academic.oup.com/mnras/article/392/1/216/1071929).

### Arguments
- `nij::Array{<:AbstractFloat,2}`: `nν`-by-`nν` symmetric noise covariance matrix, where `nν` is the number of frequency bands.
- `A::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the frequency response, for `nc` components in sky.
    - E.g., ``A = [a B]`` where `a = ones(nν)` for CMB and `B` is a `nν`-by-`nc-1` matrix for the frequency response of foreground components.
- `d::Array{<:AbstractFloat,1}`: data vector for a given pixel (or any other appropriate domain). The number of elements is `nν`.

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
- [examples/ILCPipelineSOSAT.jl](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/ILCPipelineSOSAT.jl)
  - This code applies `ilc_weights()` in harmonic domain to produce power spectra of clean polarisation maps of the CMB.
  - The code requires [Healpix.jl](https://github.com/ziotom78/Healpix.jl) and [Libsharp.jl](https://github.com/ziotom78/Libsharp.jl). It also calls python wrappers [pymaster](https://github.com/LSSTDESC/NaMaster) and [classy](https://github.com/lesgourg/class_public/wiki/Python-wrapper) via `PyCall`.
  - This is a simulation pipeline for the Small Aperture Telescope (SAT) of the [Simons Observatory](https://simonsobservatory.org). The code performs the following operations:
    1. Read in and smooth foreground maps by beams of the telescope.
    2. Generate theoretical scalar- and tensor-mode power spectra of the CMB using [CLASS](https://github.com/lesgourg/class_public).
    3. Generate Gaussian random realisations of the CMB and noise.
    4. Calculate covariance matrices `cij` of the CMB, noise, and foreground.
    5. Calculate ILC weights in harmonic domain using `ilc_weights()`.
    6. Calculate power spectra of the clean CMB maps using `ilc_clean_cij()`.
    7. Show results for visual inspection, if `showresults = true`.
    8. Calculate the tensor-to-scalar ratio and its uncertainty from the simulated realisations.
    9. Write out the results (tensor-to-scalar ratios with and without residual foreground marginalisation) to a file `ilc_results_sosat.csv` and create a PDF figure `clbb_sim_sosat.pdf` showing the cleaned B-mode power spectrum, minus noisebias, and minus the residual foreground.
  - For your reference, the results from 300 realisations are given in [results/ilc_results_sosat_300sims.csv](https://github.com/komatsu5147/CleanCMB.jl/tree/master/examples/results/ilc_results_sosat_300sims.csv). You can compute the mean and standard deviation of the tensor-to-scalar ratios and compare with the results given in Table 4 ("ILC" column) of Simons Observatory [forecast paper](https://arxiv.org/abs/1808.07445).
    - Without FG marginalisation: r = (2.0 ± 1.7) x 10<sup>-3</sup>. With marginalisation: r = (-0.9 ± 2.9) x 10<sup>-3</sup>.
