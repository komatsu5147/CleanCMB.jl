# CleanCMB

This package contains functions to enable extraction of clean maps of the cosmic microwave background (CMB). Some functions can be used to extract non-CMB astrophysical components as well.

Different algorithms exist for extraction of clean maps of the CMB (as well as of astrophysical components). The package currently supports the following algorithms:

## Internal Linear Combination (ILC)
- `ilc_weights(cij[, e, ℓid=3])`: return weights for the internal linear combination (ILC) method, following Equation (12) of [Tegmark et al., PRD, 68, 123523 (2003)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.68.123523).
- `cilc_weights(cij, a, b[, ℓid=3])`: return weights for the constrained internal linear combination (CILC) method, following Equation (20) of [Remazeilles, et al., MNRAS, 410, 2481 (2011)](https://academic.oup.com/mnras/article/410/4/2481/1007333).
- `milca_weights(cij, a, B[, ℓid=3])`: return weights for the modified internal linear combination algorithm (MILCA) method, following Equation (15) of [Hurier, et al., A&A, 558, A118 (2013)](https://www.aanda.org/articles/aa/abs/2013/10/aa21891-13/aa21891-13.html).
- Note: These papers define weights in various domain, including harmonic, wavelet, and pixel domain. You can choose to work in any domain by supplying the covariance matrix `cij` in the appropriate domain.
- `ilc_clean_cij(cij, w)`: return a covariance matrix multiplied by weights, `w' * cij * w`, with `w` returned by any of the above functions, e.g., `w = cilc_weights(cij, a, b)`. This would be a variance of the extracted component if `cij` were the same as that used in `ilc_weights()`, `cilc_weights()` or `milca_weights()`. On the other hand if, for example, `cij`, is a noise covariance matrix, this function returns the noise variance of the cleaned map for each multipole, band-power bin, pixel, etc.

### Arguments
- `cij::Array{<:AbstractFloat,2}`: `nν`-by-`nν` symmetric covariance matrix, where `nν` is the number of frequency bands.
  - or, `cij::Array{<:AbstractFloat,3}`:symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default). Here, `nℓ` is the number of elements in the relevant domain, e.g., multipoles, band-power bins, pixels, etc. `ilc_weights()`, `cilc_weights()` and `milca_weights()` are multiple-dispatch functions, which detect the dimension of the input `cij` automatically.
- `a::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. E.g., `a=[1,...,1]` for CMB. The number of elements is `nν`.
- `b::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be nulled. The number of elements is `nν`.
- `B::Array{<:AbstractFloat,2}`: `nν`-by-`nrc` matrix of the frequency response, for `nrc` components to be nulled.

### Optional arguments
- `e::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. The default value is `e=[1,...,1]`, i.e., CMB. The number of elements is `nν`.
- `ℓid::Integer=3`: location of the index for the `nℓ` domain, if the dimension of the input `cij` is `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)`.
`ℓid=1` if `cij[nℓ,nν,nν]`, `ℓid=2` if `cij[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cij[nν,nν,nℓ]`.
