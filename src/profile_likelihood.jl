"""
    loglike_beta(nij, A, d)

Profile likelihood for the frequency response vectors.

*Reference*: Equation (9) of Stompor et al., MNRAS, 392, 216 (2009)

# Arguments
- `nij::Array{<:AbstractFloat,2}`: symmetric noise covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `A::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the frequency response, for `nc` components in sky.
    - E.g., ``A = [a B]`` where `a = ones(nν)` for CMB and `B` is a `nν`-by-`nc-1` matrix for the frequency response of foreground components.
- `d::Array{<:AbstractFloat,1}`: data vector for a given pixel (or any other appropriate domain). The number of elements is `nν`.
"""
function loglike_beta(
    nij::Array{T,2},
    A::Array{T,2},
    d::Array{T,1},
) where {T<:AbstractFloat}
    if size(nij)[1] ≠ size(nij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    else
        M = Symmetric(nij)
    end
    x = M \ d
    s = inv(A' * inv(M) * A) * A'x
    lnlike = 0.5 * x'A * s
end

"""
    loglike_beta_deriv(nij, A, dAdβ, d)

Derivative of the profile likelihood for the frequency response vectors with respect to a foreground parameter

*Reference*: Equation (A1) of Errard et al., PRD, 84, 063005 (2011)

# Arguments
- `nij::Array{<:AbstractFloat,2}`: symmetric noise covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `A::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the frequency response, for `nc` components in sky.
    - E.g., ``A = [a B]`` where `a = ones(nν)` for CMB and `B` is a `nν`-by-`nc-1` matrix for the frequency response of foreground components.
- `dAdβ::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the derivative of the frequency response with respect to a foreground parameter.
    - E.g., ``dAdβ = [zeros(nν) dsynch/dβs zeros(nν)]`` where `zeros(nν)` for CMB and dust because they do not depend on the synchrotron index `βs`.
- `d::Array{<:AbstractFloat,1}`: data vector for a given pixel (or any other appropriate domain). The number of elements is `nν`.
"""
function loglike_beta_deriv(
    nij::Array{T,2},
    A::Array{T,2},
    dAdβ::Array{T,2},
    d::Array{T,1},
) where {T<:AbstractFloat}
    if size(nij)[1] ≠ size(nij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    else
        M = Symmetric(nij)
    end
    x = M \ d
    s = inv(A' * inv(M) * A) * A'x
    lnlike_deriv = s' * dAdβ' * (x - M \ (A * s))
end
