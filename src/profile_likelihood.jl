"""
    loglike_beta(nij, A, cij)

Profile likelihood for the frequency response vectors.

*Reference*: Based on Equation (9) of Stompor et al., MNRAS, 392, 216 (2009)

# Arguments
- `nij::Array{<:AbstractFloat,2}`: symmetric noise covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `A::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the frequency response, for `nc` components in sky.
    - E.g., ``A = [a B]`` where `a = ones(nν)` for CMB and `B` is a `nν`-by-`nc-1` matrix for the frequency response of foreground components.
- `cij::Array{<:AbstractFloat,2}`: symmetric covariance matrix with the dimention of `(nν, nν)`.
"""
function loglike_beta(
    nij::Array{T,2},
    A::Array{T,2},
    cij::Array{T,2},
) where {T<:AbstractFloat}
    if size(nij)[1] ≠ size(nij)[2] || size(cij)[1] ≠ size(cij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    elseif size(nij)[1] ≠ size(A)[1] || size(nij)[1] ≠ size(cij)[1]
        throw(DimensionMismatch("dimensions of input matrices do not match"))
    else
        N = Symmetric(nij)
        M = Symmetric(cij)
    end
    Ninv = inv(N)
    lnlike = 0.5 * tr(Ninv * A * inv(A' * Ninv * A) * A' * Ninv * M)
end

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
    elseif size(nij)[1] ≠ size(A)[1] || size(nij)[1] ≠ length(d)
        throw(DimensionMismatch("dimensions of matrices/vector do not match"))
    else
        M = Symmetric(nij)
    end
    x = M \ d
    s = inv(A' * inv(M) * A) * A'x
    lnlike = 0.5 * x'A * s
end

"""
    loglike_beta_grad(nij, A, dAdβ, d)

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
function loglike_beta_grad(
    nij::Array{T,2},
    A::Array{T,2},
    dAdβ::Array{T,2},
    d::Array{T,1},
) where {T<:AbstractFloat}
    if size(nij)[1] ≠ size(nij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    elseif size(nij)[1] ≠ size(A)[1] ||
           size(nij)[1] ≠ length(d) ||
           size(A)[1] ≠ size(dAdβ)[1] ||
           size(A)[2] ≠ size(dAdβ)[2]
        throw(DimensionMismatch("dimensions of matrices/vector do not match"))
    else
        M = Symmetric(nij)
    end
    x = M \ d
    s = inv(A' * inv(M) * A) * A'x
    lnlike_grad = s' * dAdβ' * (x - M \ (A * s))
end

"""
    loglike_beta_hess(nij, A, dAdβI, dAdβJ, d)

Hessian of the profile likelihood for the frequency response vectors with respect to foreground parameters

*Reference*: Equation (5) of Errard et al., PRD, 84, 063005 (2011)

# Arguments
- `nij::Array{<:AbstractFloat,2}`: symmetric noise covariance matrix with the dimention of `(nν, nν)` where `nν` is the number of frequency bands.
- `A::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the frequency response, for `nc` components in sky.
    - E.g., ``A = [a B]`` where `a = ones(nν)` for CMB and `B` is a `nν`-by-`nc-1` matrix for the frequency response of foreground components.
- `dAdβI::Array{<:AbstractFloat,2}` and `dAdβJ::Array{<:AbstractFloat,2}`: `nν`-by-`nc` matrix of the derivative of the frequency response with respect to a foreground parameter.
    - E.g., ``dAdβI = [zeros(nν) dsynch/dβs zeros(nν)]`` where `zeros(nν)` for CMB and dust because they do not depend on the synchrotron index `βs`.
- `d::Array{<:AbstractFloat,1}`: data vector for a given pixel (or any other appropriate domain). The number of elements is `nν`.
"""
function loglike_beta_hess(
    nij::Array{T,2},
    A::Array{T,2},
    dAdβI::Array{T,2},
    dAdβJ::Array{T,2},
    d::Array{T,1},
) where {T<:AbstractFloat}
    if size(nij)[1] ≠ size(nij)[2]
        throw(DimensionMismatch("covariance matrix must be a square matrix"))
    elseif size(nij)[1] ≠ size(A)[1] ||
           size(nij)[1] ≠ length(d) ||
           size(A)[1] ≠ size(dAdβI)[1] ||
           size(A)[2] ≠ size(dAdβI)[2] ||
           size(A)[1] ≠ size(dAdβJ)[1] ||
           size(A)[2] ≠ size(dAdβJ)[2]
        throw(DimensionMismatch("dimensions of matrices/vector do not match"))
    else
        M = Symmetric(nij)
    end
    x = M \ d
    AtMinvAinv = inv(A' * inv(M) * A)
    s = AtMinvAinv * A'x
    u = M \ (dAdβJ * s)
    v = M \ (A * AtMinvAinv * A'u)
    lnlike_hess = tr(dAdβI' * (v - u) * s')
end
