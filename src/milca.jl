"""
    milca_weights(cij, a, B)

This function returns weights (a vector of the number of frequency channels) of the Modified ILC Algorith (MILCA) method.

*Reference*: Equation (15) of Hurier, Macías-Pérez, Hildebrandt, A&A, 558, A118 (2013)

# Arguments
- `cij::Array{<:AbstractFloat,2}`: symmetric covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `a::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. E.g., `a = [1,...,1]` for CMB.
- `B::Array{<:AbstractFloat,2}`: `nν`-by-`nrc` matrix of the frequency response, for `nrc` components to be nulled.
"""
function milca_weights(
    cij::Array{T,2},
    a::Array{T,1},
    B::Array{T,2},
) where {T<:AbstractFloat}
    if size(cij)[1] ≠ size(cij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    elseif size(cij)[1] ≠ length(a) || size(cij)[1] ≠ size(B)[1]
        throw(DimensionMismatch("dimensions of the covariance matrix and the frequency response vector or matrix do not match"))
    else
        M = Symmetric(cij)
    end
    e = zeros(size(B)[2] + 1)
    e[1] = 1
    F = [a B]
    w = inv(M) * F * inv(F' * inv(M) * F) * e
end

"""
    milca_weights(cijℓ, a, B[, ℓid=3])

This function returns weights (a `nν`-by-`nℓ` matrix) of the Modified ILC Algorith (MILCA) method.

Here, `nν` is the number of frequency channels and `nℓ` is the number of elements in the relevant domain, e.g., multipoles, band-power bins, pixels, etc.

*Reference*: Equation (15) of Hurier, Macías-Pérez, Hildebrandt, A&A, 558, A118 (2013)

# Arguments
- `cijℓ::Array{<:AbstractFloat,3}`: symmetric covariance matrix with the dimention of `(nℓ, nν, nν)`, `(nν, nℓ, nν)` or `(nν, nν, nℓ)` (default).
- `a::Array{<:AbstractFloat,1}`: vector of the frequency response, for the component to be extracted. E.g., `a = [1,...,1]` for CMB.
- `B::Array{<:AbstractFloat,2}`: `nν`-by-`nrc` matrix of the frequency response, for `nrc` components to be nulled.

# Optional Arguments
- `ℓid::Integer=3`: location of the index for the `nℓ` domain. `ℓid=1` if `cijℓ[nℓ,nν,nν]`, `ℓid=2` if `cijℓ[nν,nℓ,nν]`, and `ℓid=3` (the default value) if `cijℓ[nν,nν,nℓ]`.
"""
function milca_weights(
    cijℓ::Array{T,3},
    a::Array{T,1},
    B::Array{T,2},
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
    if length(a) ≠ nν || size(B)[1] ≠ nν
        throw(DimensionMismatch("dimensions of the covariance matrix and the frequency response vector/matrix do not match"))
    end
    wℓ = zeros(nν, nℓ) # ILC weights
    e = zeros(size(B)[2] + 1)
    e[1] = 1
    F = [a B]
    for iℓ = 1:nℓ
        if ℓid == 3
            M = Symmetric(cijℓ[:, :, iℓ])
        elseif ℓid == 2
            M = Symmetric(cijℓ[:, iℓ, :])
        else
            M = Symmetric(cijℓ[iℓ, :, :])
        end
        wℓ[:, iℓ] = inv(M) * F * inv(F' * inv(M) * F) * e
    end
    return wℓ
end
