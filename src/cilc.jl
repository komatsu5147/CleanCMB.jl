"""
    cilc_weights(cij, a, b)

This function returns weights (a vector of the number of frequency channels) of the constrained ILC (CILC) method.

*Reference*: Equation (20) of Remazeilles, Delabrouille, Cardoso, MNRAS, 410, 2481 (2011)

# Arguments
- `cij::Array{<:AbstractFloat,2}`: symmetric covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `a::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. E.g., a=[1,...,1] for CMB.
- `b::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be nulled.
"""
function cilc_weights(
    cij::Array{T,2},
    a::Array{T,1},
    b::Array{T,1},
) where {T<:AbstractFloat}
    if size(cij)[1] ≠ size(cij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    elseif size(cij)[1] ≠ length(a) || size(cij)[1] ≠ length(b)
        throw(DimensionMismatch("dimensions of the covariance matrix and the frequency response vector do not match"))
    elseif a ≈ b
        throw(ErrorException("vectors of the frequency response are too similar"))
    else
        M = Symmetric(cij)
    end
    x = M \ a # R^-1 a
    y = M \ b # R^-1 b
    w = (b'y * x - a'y * y) / (a'x * b'y - (a'y)^2) # ILC weights
end

"""
    cilc_weights(cijℓ, a, b[, ℓid=3])

This function returns weights (a `nν`-by-`nℓ` matrix) of the constrained ILC (CILC) method.

Here, `nν` is the number of frequency channels and `nℓ` is the number of elements in the relevant domain, e.g., multipoles, band-power bins, pixels, etc.

*Reference*: Equation (20) of Remazeilles, Delabrouille, Cardoso, MNRAS, 410, 2481 (2011)

# Arguments
- `cijℓ::Array{<:AbstractFloat,3}`: symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default).
- `a::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. E.g., a=[1,...,1] for CMB.
- `b::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be nulled.

# Optional Arguments
- `ℓid::Integer=3`: location of the index for the `nℓ` domain. `ℓid=1` if `cijℓ[nℓ,nν,nν]`, `ℓid=2` if `cijℓ[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cijℓ[nν,nν,nℓ]`.
"""
function cilc_weights(
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
    elseif a ≈ b
        throw(ErrorException("vectors of the frequency response are too similar"))
    else
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
    cilc_clean_cij(cij, w)

This function returns power of the extracted component for a given element of the relevant domain, e.g., multipole, band-power bin, pixel, etc.

*Reference*: Tegmark et al., Phys. Rev. D, 68, 123523 (2003)

# Arguments
- `cij::Array{<:AbstractFloat,2}`: symmetric covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `w::Array{<:AbstractFloat,1}`: ILC weights.
"""
function cilc_clean_cij(cij::Array{T,2}, w::Array{T,1}) where {T<:AbstractFloat}
    cl = ilc_clean_cij(cij, w)
end

"""
    cilc_clean_cij(cijℓ, wℓ[, ℓid::Integer=3])

This function returns a vector of the power of the extracted component, with elements in the relevant domain, e.g., multipole, band-power bin, pixel, etc.

*Reference*: Tegmark et al., Phys. Rev. D, 68, 123523 (2003)

# Arguments
- `cijℓ::Array{<:AbstractFloat,3}`: symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default).
- `wℓ::Array{<:AbstractFloat,2}`: ILC weights (a `nν`-by-`nℓ` matrix).
    - Here, `nν` is the number of frequency bands and `nℓ` is the number of elements in the relevant domain.

# Optional Arguments
- `ℓid::Integer=3`: location of the index for the `nℓ` domain. `ℓid=1` if `cijℓ[nℓ,nν,nν]`, `ℓid=2` if `cijℓ[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cijℓ[nν,nν,nℓ]`.
"""
function cilc_clean_cij(
    cijℓ::Array{T,3},
    wℓ::Array{T,2},
    ℓid::Integer = 3,
) where {T<:AbstractFloat}
    cl = ilc_clean_cij(cijℓ, wℓ, ℓid)
end
