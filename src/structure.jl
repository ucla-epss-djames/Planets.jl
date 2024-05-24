using SpecialFunctions: gamma
using PhysicalConstants.CODATA2014: G
using QuadGK: quadgk

export dmdr, dPdr
export planet_m, planet_g, planet_mmotion, planet_mu, planet_cmu

"""
    dmdr(r::Float64, ρ::Float64)::Float64

Mass in a spherical shell.

# Arguments
- `r::Float64` - radius
- `ρ::Float64` - density
"""
dmdr(r::Float64, ρ::Float64)::Float64 = 4π * ρ * r^2

"""
    dPdr(r::Float64, ρ::Float64, gm::Float64)::Float64

Pressure due to hydrostatic equilibrium. The gravitational mass is Newton's
Gravitational constant multiplied by the mass of the body.

# Arguments
- `r::Float64`  - radius
- `ρ::Float64`  - density
- `gm::Float64` - gravitational mass
"""
dPdr(r::Float64, ρ::Float64, gm::Float64)::Float64 = -gm * ρ / r^2

"dPdr(ρ::Float64, g::Float64)::Float64"
dPdr(ρ::Float64, g::Float64)::Float64 = -ρ * g

"""
    planet_g(r::Float64, gm::Float64)::Float64

Gravity due to a spherical mass.

# Arguments
- `r::Float64`  - radius
- `gm::Float64` - gravitational mass
"""
planet_g(r::Float64, gm::Float64)::Float64 = r == 0 ? 0 : gm / r^2

"""
    calc_gravity(r0::Float64, r1::Float64, m::Float64,
                      ρ::Function)::NTuple{2, Float64}

Calculates the gravity of a planetary body.

# Arguments
- `r0::Float64` - initial radius
- `r1::Float64` - final radius
- `m::Float64` - mass
- `ρ::Function` - density function
"""
function calc_gravity(r0::Float64, r1::Float64, m::Float64,
                      ρ::Function)::NTuple{2, Float64}

    res, err = quadgk(x -> dmdr(x, ρ(x)), r0, r1)
    m += res

    g = planet_g(r1, m*G.val)

    return (m, g)
end

"""
    calc_gravity(r0::Float64, r1::Float64, m::Float64,
                      ρ::Float64)::Ntuple{2, Float64}

Calculates the gravity of a planetary body.

# Arguments
- `r0::Float64` - initial radius
- `r1::Float64` - final radius
- `m::Float64` - mass
- `ρ::Float64` - density
"""
function calc_gravity(r0::Float64, r1::Float64, m::Float64,
                      ρ::Float64)::Ntuple{2, Float64}

    res, err = quadgk(x -> dmdr(x, ρ), r0, r1)
    m += res

    g = planet_g(r1, m*G.val)

    return (m, g)
end


"""
    planet_mmotion(gm::Float64, a::Float64)::Float64

Mean motion due to a planetary body.

# Arguments
- `gm::Float64` - gravitational mass of planet
- `a::Float64`  - radial orbit of moon
"""
planet_mmotion(gm::Float64, a::Float64)::Float64 = a == 0 ? 0 : sqrt(gm / a^3)

"""
    planet_mu(r::Float64, g::Float64, ρ::Float64)::Float64

Veritcal stress due to gravity.

# Arguments
- `r::Float64` - radius
- `g::Float64` - gravity
- `ρ::Float64` - density
"""
planet_mu(r::Float64, g::Float64, ρ::Float64)::Float64 = ρ * g * r


"""
    cmu_maxwell(μ::Float64, ω::Float64, η::Float64)::ComplexF64

Maxwell rheaology using the tidal forcing frequency. Refer to Turcotte &
Schubert 2002 or Storch & Lai 2014 for a decription of this rheaology.

# Arguments
- `μ::Float64` - shear modulus (rigidity)
- `ω::Float64` - rotational frequency
- `η::Float64` - viscosity
"""
function cmu_maxwell(μ::Float64, ω::Float64, η::Float64)::ComplexF64

    ω_m = μ / η

    cmu = μ / (1 - (ω_m/ω)*im)

    return cmu

end

"""
    cmu_SLS(μ0::Float64, ω::Float64, η::Float64,
            μ_f::Union{Int64, Float64})::ComplexF64

Standard Linear Solid model. Refer to Stixrude et al. 2021 or Norwick & Berry
1972 for a description of this rheaology.

# Arguments
- `μ0::Float64`  - unrelaxed shear modulus
- `ω::Float64`   - rotational frequency
- `η::Float64`   - viscosity
- `μ_f::Float64` - SLS parameter s.t. μ₁/μ₀
"""
function cmu_SLS(μ0::Float64, ω::Float64, η::Float64,
                 μ_f::Union{Int64, Float64})::ComplexF64

    μ1 = μ0 * μ_f

    dμ = μ0 * (1 - μ1 / (μ0 + μ1))

    τ = η / (μ0 + μ1)

    cmu = μ0 - dμ / (1 + τ*ω*im)

end

"""
    cmu_andrade(μ::Float64, ω::Float64, η::Float64, α::Float64;
                β::Float64=0.)::ComplexF64

Andrade model. Refer to Castillo-Rogez et al. 2011 or Dumoulin et al. 2017 for a
description of this rheaology.

# Arguments
- `μ::Float64` - shear modulus
- `ω::Float64` - rotational frequency
- `η::Float64` - viscosity
- `α::Float64` - frequency dependence (Andrade exponent)
- `β::Float64` - amplitude of the transient response
"""
function cmu_andrade(μ::Float64, ω::Float64, η::Float64, α::Float64;
                     β::Float64=0.)::ComplexF64

    if(β == 0) β = μ^(α - 1) * η^-α end

    J = 1 / μ - im / (η * ω) + β*(ω*im)^-α * gamma(1 + α)

    cmu = 1 / J

end

"""
    planet_cmu(μ::Float64, ω::Float64, η::Float64, r::Float64, g::Float64,
               ρ::Float64, μ_f::NTuple{3, Union{Int64, Float64}},
               model::Int64)::ComplexF64

Calculates the complex shear modulus (CMU) based on model input desire.
Viscosity will determine if the CMU model is used or if the infinite limit is
used. If viscosity is 0.0, a small shear modulus is produced for stability.
Refer to `cmu_maxwell`, `cmu_SLS`, or `cmu_andrade` for infromation on the
models.

# Arguments
- `μ::Float64` - shear modulus
- `ω::Float64` - rotational frequency
- `η::Float64` - viscosity
- `r::Float64` - radius
- `g::Float64` - gravity
- `ρ::Float64` - density
- `μ_f::NTuple{3, Union{Int64, Float64}}` - rheaology parameters
                                            ([1] SLS factor, [2,3] Andrade alpha, beta)
- `model::Int64` - model desire ([1] Maxwell, [2] SLS, [3] Andrade)
"""
function planet_cmu(μ::Float64, ω::Float64, η::Float64, r::Float64, g::Float64,
                    ρ::Float64, μ_f::NTuple{3, Union{Int64, Float64}},
                    model::Int64)::ComplexF64

    if η == 0.0
        # if eta is small, produce small shear

        smu = 1e-4 * planet_mu(r, g, ρ)
        cmu = smu + 0*im

    elseif η == -1
        # if eta is infinite, produce the real shear

        cmu = μ + 0*im

    else
        # else use the desired rheaology

        if model == 1
            cmu = cmu_maxwell(μ, ω, η)
        elseif model == 2
            cmu = cmu_SLS(μ, ω, η, μ_f[1])
        elseif model == 3
            cmu = cmu_andrade(μ, ω, η, μ_f[2], β=μ_f[3])
        end

    end

    if(real(cmu) == 0.0)
        # final check if shear is 0.0

        cmu = 1e-9 + 0*im

    end

    return cmu

end
