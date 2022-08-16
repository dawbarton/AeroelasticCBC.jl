export unsteady_model!, unsteady_model, unsteady_model_case1, unsteady_model_case2

const BASE_PARS = ComponentArray(; U = 15, b = 0.15, a = -0.5, ρ = 1.204, m_w = 5.3,
                                 m_T = 16.9)

"""
    unsteady_model_case1()

Return a vector (`ComponentVector`) of model parameters for case 1.
"""
function unsteady_model_case1()
    return ComponentArray(BASE_PARS;
                          I_α = 0.1724,
                          c_α = 0.5628,
                          c_h = 14.5756,
                          k_α = 54.1162,
                          k_α₂ = 751.6,
                          k_α₃ = 5006.7,
                          k_h = 3.5294e+03,
                          x_α = 0.24)
end

"""
    unsteady_model_case2()

Return a vector (`ComponentVector`) of model parameters for case 2.
"""
function unsteady_model_case2()
    return ComponentArray(BASE_PARS;
                          I_α = 0.1726,
                          c_α = 1.0338,
                          c_h = 15.4430,
                          k_α = 60.291,
                          k_α₂ = 774.7,
                          k_α₃ = 3490.7,
                          k_h = 3.3178e+03,
                          x_α = 0.234)
end

# using Symbolics
# @variables b a ρ m_w m_T I_α c_α c_h k_α k_α₂ k_α₃ k_h x_α U π
# M = [
#     m_T+π*ρ*b^2         m_w*x_α*b-a*π*ρ*b^3   0
#     m_w*x_α*b-a*π*ρ*b^3 I_α+π*(1/8+a^2)*ρ*b^4 0
#     0                   0                     1
# ]
# M⁻¹ = inv(M)
# repr(M⁻¹)

"""
    unsteady_model!(du::ComponentArray, u::ComponentArray, p, t)

Return (in `du`) the right-hand side of the differential equations corresponding to the
unsteady aerodynamic model. `du`, `u` and `p` are accessed by properties (e.g., `u.α`) and
so ordering within the input/output vectors is arbitrary.
"""
function unsteady_model!(du::ComponentArray, u::ComponentArray, p, t)
    # Unpack state variables (as required)
    v = @SVector [u.h, u.α, u.w]
    v′ = @SVector [u.h′, u.α′, u.w′]

    # Unpack parameters
    (; U, b, a, ρ, m_w, m_T, I_α, c_α, c_h, k_α, k_α₂, k_α₃, k_h, x_α) = p  # order doesn't matter since this gets expanded to U = p.U, etc

    # Sears coefficients (Abdelkefi et al, Nonlinear Dynamics 2013)
    c₀ = 1
    c₁ = 0.1650
    c₂ = 0.0455
    c₃ = 0.335
    c₄ = 0.3

    # Build linear terms
    M⁻¹ = @SMatrix [(true + (((b * m_w * x_α - a * π * ρ * (b^3)) / (m_T + π * ρ * (b^2))) * (b * m_w * x_α - a * π * ρ * (b^3))) / (I_α + (-((b * m_w * x_α - a * π * ρ * (b^3))^2)) / (m_T + π * ρ * (b^2)) + π * ρ * (b^4) * (0.125 + a^2)))/(m_T + π * ρ * (b^2)) (-((b * m_w * x_α - a * π * ρ * (b^3)) / (I_α + (-((b * m_w * x_α - a * π * ρ * (b^3))^2)) / (m_T + π * ρ * (b^2)) + π * ρ * (b^4) * (0.125 + a^2))))/(m_T + π * ρ * (b^2)) 0
                    (-((b * m_w * x_α - a * π * ρ * (b^3)) / (m_T + π * ρ * (b^2))))/(I_α + (-((b * m_w * x_α - a * π * ρ * (b^3))^2)) / (m_T + π * ρ * (b^2)) + π * ρ * (b^4) * (0.125 + a^2)) true/(I_α + (-((b * m_w * x_α - a * π * ρ * (b^3))^2)) / (m_T + π * ρ * (b^2)) + π * ρ * (b^4) * (0.125 + a^2)) 0
                    0.0 0.0 1.0]

    D = @SMatrix [c_h+2 * π * ρ * b * U * (c₀ - c₁ - c₃) (1+(c₀ - c₁ - c₃) * (1 - 2a))*π*ρ*b^2*U 2π*ρ*U^2*b*(c₁ * c₂+c₃ * c₄)
                  -2π*(a+1 / 2)*ρ*b^2*(c₀ - c₁-c₃)*U c_α+(1 / 2 - a) * (1 - (c₀ - c₁ - c₃) * (1 + 2a)) * π * ρ * b^3 * U -2π*ρ*b^2*U^2*(a+1 / 2)*(c₁ * c₂+c₃ * c₄)
                  -1/b a-1 / 2 (c₂ + c₄) * U/b]

    K = @SMatrix [k_h 2π*ρ*b*U^2*(c₀ - c₁-c₃) 2π*ρ*U^3*c₂*c₄*(c₁+c₃)
                  0 k_α-2π * (1 / 2 + a) * ρ * (c₀ - c₁ - c₃) * b^2 * U^2 -2π*ρ*b*U^3*(a+1 / 2)*c₂*c₄*(c₁+c₃)
                  0 -U/b c₂ * c₄ * U^2/b^2]

    # Build nonlinear terms
    α = u.α
    N = @SVector [0, k_α₂ * α^2 + k_α₃ * α^3, 0]

    # Construct equations
    v″ = -M⁻¹ * (D * v′ + K * v + N)

    # Reassemble
    du.h = v′[1]
    du.α = v′[2]
    du.w = v′[3]
    du.h′ = v″[1]
    du.α′ = v″[2]
    du.w′ = v″[3]

    return du
end

"""
    unsteady_model!(du, u, p, t)

Return (in `du`) the right-hand side of the differential equations corresponding to the
unsteady aerodynamic model. The ordering of the state vector is assumed to be
`(h, h′, α, α′, w, w′)` and is internally converted to a `ComponentArray`.
"""
function unsteady_model!(du, u, p, t)
    ax = ComponentArrays.Axis((h = 1, h′ = 2, α = 3, α′ = 4, w = 5, w′ = 6))
    unsteady_model!(ComponentArray(du, ax), ComponentArray(u, ax), p, t)
    return du
end

"""
    unsteady_model(u, p)

A wrapper for `unsteady_model!` for use with `BifurcationKit.jl`.
"""
unsteady_model(u, p) = unsteady_model!(similar(u), u, p, 0)
