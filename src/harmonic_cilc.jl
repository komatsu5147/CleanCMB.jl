"""
    harmonic_cilc_weights(cij, a, b)

This function returns weights (a vector of the number of frequency channels) of the harmonic constrained ILC (CILC) method for a given multipole or band-power bin.

*Reference*: Equation (20) of Remazeilles, Delabrouille, Cardoso, MNRAS, 410, 2481 (2011)

# Arguments
- `cij::Array{T,2}`: symmetric covariance matrix for a given multipole (or a band-power bin) with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `a::Array{T,1}`: vector of the frequency response, which is extracted. E.g., a=[1,...,1] for CMB.
- `b::Array{T,1}`: vector of the frequency response, which is nulled.
"""
function harmonic_cilc_weights(
    cij::Array{T,2},
    a::Array{T,1},
    b::Array{T,1},
) where {T<:AbstractFloat}
    if size(cij)[1] ≠ size(cij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    elseif size(cij)[1] ≠ length(a) || size(cij)[1] ≠ length(b)
        throw(DimensionMismatch("dimensions of the covariance matrix and the frequency response vector do not match"))
    else
        M = Symmetric(cij)
    end
    x = M \ a # R^-1 a
    y = M \ b # R^-1 b
    w = (b'y * x - a'y * y) / (a'x * b'y - (a'y)^2) # ILC weights
end

"""
    harmonic_cilc_weights(cijℓ, a, b[, ℓid=3])

This function returns weights (a nν-by-nℓ matrix) of the harmonic constrained ILC (CILC) method. The first element is a frequency channel and the second element is a multipole or band-power bin.

*Reference*: Equation (20) of Remazeilles, Delabrouille, Cardoso, MNRAS, 410, 2481 (2011)

# Arguments
- `cijℓ::Array{T,3}`: symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default) where `nν` is the number of frequency bands and `nℓ` is the number of multipoles or band-power bins.
- `a::Array{T,1}`: vector of the frequency response, which is extracted. E.g., a=[1,...,1] for CMB.
- `b::Array{T,1}`: vector of the frequency response, which is nulled.

# Optional Arguments
- `ℓid::Integer=3`: location of the index for multipoles or band-power bins. `ℓid=1` if `cijℓ[nℓ,nfreq,nfreq]`, `ℓid=2` if `cijℓ[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cijℓ[nν,nν,nℓ]`.
"""
function harmonic_cilc_weights(
    cijℓ::Array{T,3},
    a::Array{T,1},
    b::Array{T,1},
    ℓid::Integer = 3,
) where {T<:AbstractFloat}
    if ℓid > 3 || ℓid < 1
        throw(DomainError(ℓid, "ℓid must be 1, 2, or 3"))
    end
    if (ℓid == 3 && size(cijℓ)[1] ≠ size(cijℓ)[2]) ||
       (ℓid == 2 && size(cijℓ)[1] ≠ size(cijℓ)[3]) ||
       (ℓid == 1 && size(cijℓ)[2] ≠ size(cijℓ)[3])
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    end
    nℓ = size(cijℓ)[ℓid]
    nν = ifelse(ℓid == 3, size(cijℓ)[1], size(cijℓ)[3])
    if length(a) ≠ nν || length(b) ≠ nν
        throw(DimensionMismatch("dimensions of the covariance matrix and the frequency response vector do not match"))
    end
    wℓ = zeros(nν, nℓ) # ILC weights
    for iℓ = 1:nℓ
        if ℓid == 3
            M = Symmetric(cijℓ[:, :, iℓ])
        elseif ℓid == 2
            M = Symmetric(cijℓ[:, iℓ, :])
        else
            M = Symmetric(cijℓ[iℓ, :, :])
        end
        x = M \ a # R^-1 a
        y = M \ b # R^-1 b
        wℓ[:, iℓ] = (b'y * x - a'y * y) / (a'x * b'y - (a'y)^2)
    end
    return wℓ
end

"""
    harmonic_cilc_clean_cl(cij, w)

This function returns a cleaned CMB power spectrum at a given multipole or a band-power bin.

*Reference*: Tegmark et al., Phys. Rev. D, 68, 123523 (2003)

# Arguments
- `cij::Array{T,2}`: symmetric covariance matrix for a given multipole (or a band-power bin) with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `w::Array{T,1}`: ILC weights for a given multipole (or a band-power bin)
"""
function harmonic_cilc_clean_cl(
    cij::Array{T,2},
    w::Array{T,1},
) where {T<:AbstractFloat}
    cl = harmonic_ilc_clean_cl(cij, w)
end

"""
    harmonic_cilc_clean_cl(cijℓ, wℓ[, ℓid::Integer=3])

This function returns a cleaned CMB power spectrum (a vector of the number of multipoles or band-power bins).

*Reference*: Tegmark et al., Phys. Rev. D, 68, 123523 (2003)

# Arguments
- `cijℓ::Array{T,3}`: symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default) where `nν` is the number of frequency bands and `nℓ` is the number of multipoles or band-power bins.
- `wℓ::Array{T,2}`: ILC weights (a nν-by-nℓ matrix)

# Optional Arguments
- `ℓid::Integer=3`: location of the index for multipoles or band-power bins. `ℓid=1` if `cijℓ[nℓ,nfreq,nfreq]`, `ℓid=2` if `cijℓ[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cijℓ[nν,nν,nℓ]`.
"""
function harmonic_cilc_clean_cl(
    cijℓ::Array{T,3},
    wℓ::Array{T,2},
    ℓid::Integer = 3,
) where {T<:AbstractFloat}
    cl = harmonic_ilc_clean_cl(cijℓ, wℓ, ℓid)
end
