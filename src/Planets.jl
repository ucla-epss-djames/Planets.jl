# structure based off of Stixrude '21 et al.
# design still being added to it
# NOTE:
# - add mutator for mass and gm --> may be in PlanetEvolution.jl
# - update matrix containers --> Tidal.jl
# - is Density.jl still a thing?
# - move BM3 and other extrapolation schemes --> Elastic.jl
# - update packages across all these modules
# - push all repos by the end of the quarter
# - also update BlockingMethod with changes Lars made --> BlockingMethod.jl
# - Moon, retro flag can perhaps be removed since structure code will be added
# - look back on Andrade model
module Planets

include("structure.jl")
include("thermal.jl")

using PhysicalConstants.CODATA2014: σ

export Planet, Moon
export lumin_internal, lumin_core

"""
    Planet

Planet struct used for PlanetEvolution.jl.

# Fields
- `name       :: String`                - name of planet
- `layers     :: Int64`                 - number of layers in structure file
- `R          :: Float64`               - radius of planet
- `M          :: Float64`               - mass of planet
- `GM         :: Float64`               - gravitational mass of planet
- `ω          :: Float64`               - rotational frequency
- `C_p        :: Union{Int64, Float64}` - specific heat
- `α          :: Float64`               - thermal expansivity
- `k          :: Float64`               - thermal conductivity
- `κ          :: Float64`               - thermal diffusivity
- `T0         :: Union{Int64, Float64}` - phase transition reference T
- `P0         :: Union{Int64, Float64}` - phase transition reference P
- `a          :: Float64`               - Simon pressure
- `b          :: Float64`               - Simon gradient
- `B          :: Float64`               - effective temperature exponent
- `∇          :: Float64`               - adiabtatic gradient
- `T_eq       :: Union{Int64, Float64}` - radiative equilibrium temperature
- `T_ef       :: Union{Int64, Float64}` - Today's effective temperature
- `T1         :: Union{Int64, Float64}` - Temp at 1 bar
- `η0         :: Float64`               - reference viscosity
- `A          :: Int64`                 - viscosity exponent
- `Ra         :: Int64`                 - critical Rayleigh number
- `μ_f        :: NTuple{3, Union{Int64, Float64}}`
                               - rhealogoy parameters
                                  (1 - SLS(mu factor),
                                   2,3 - Andrade(alpha, beta))
- `rhea_model :: Int64`        - modulus model(1:maxwell, 2:SLS, 3:Andrade)
"""
struct Planet

    name :: String

    # structure parameters
    layers :: Int64
    R      :: Float64
    M      :: Float64
    GM     :: Float64

    # orbital parameters
    ω :: Float64

    # thermal parameters
    C_p  :: Union{Int64, Float64}
    α    :: Float64
    k    :: Float64
    κ    :: Float64
    T0   :: Union{Int64, Float64}
    P0   :: Union{Int64, Float64}
    a    :: Float64
    b    :: Float64
    B    :: Float64
    ∇    :: Float64
    T_eq :: Union{Int64, Float64}
    T_ef :: Union{Int64, Float64}
    T1   :: Union{Int64, Float64}

    # rheaology parameters
    η0         :: Float64
    A          :: Int64
    Ra         :: Int64
    μ_f        :: NTuple{3, Union{Int64, Float64}}
    rhea_model :: Int64

end

"""
    Moon

Moon struct used for PlanetEvolution.jl. Primarily for the tidal calculation and
orbital migration.

# Fields
- `name  :: String`  - name of moon
- `m     :: Float64` - mass of moon
- `gm    :: Float64` - gravitational mass of moon
- `a     :: Float64` - orbital radius
- `i     :: Float64` - inclination
- `retro :: Bool`    - is moon in retrograde motion
"""
struct Moon

    name :: String

    # structure parameters
    m  :: Float64
    gm :: Float64

    # orbital parameters
    a :: Float64
    i :: Float64

    # flags
    retro :: Bool

end

"""
    lumin_internal(R::Float64, T_ef::Union{Int64, Float64},
                   T_eq::Union{Int64, Float64})::Float64

Calculates the total luminosity of the interior.

# Arguments
- `R::Float64`    - radius
- `T_ef::Union{Int64, Float64}` - effective temperature
- `plnt::Planet` - Planet struct for radiative equilibrium temperature
"""
function lumin_internal(R::Float64, T_ef::Union{Int64, Float64},
                        plnt::Planet)::Float64
    4π*R^2 * σ.val * (T_ef^4 - plnt.T_eq^4)
end

"""
    lumin_core(T1::Union{Int64, Float64}, Ti::Float64, c::Float64,
               P_c::Float64, T_c::Float64, ρ_c::Float64, g_c::Float64,
               P1::Float64,
               plnt::NamedTuple{(:α, :k, :κ, :T0, :P0, :a, :b, :∇, :η0, :A, :Ra),
                                Float64})::NamedTuple{(:L_c, :K, :ΔT, :Γ, :Ra), Float64}

Calculates the luminosty due to the core. Refer to Stixrude et al. 2021 (eq 9-16).

# Arguments
- `T1::Union{Int64, Float64}`  - temperature of envelope
- `Ti::Float64`  - temperature at top of thermal boundary
- `c::Float64`   - radius of core
- `P_c::Float64` - pressure at top of the core
- `T_c::Float64` - temperature at top of the core
- `ρ_c::Float64` - density at top of the core
- `g_c::Float64` - gravity at top of the core
- `P1::Float64`  - reference pressure
- `plnt::Planet` - planet struct
"""
function lumin_core(T1::Union{Int64, Float64}, Ti::Float64, c::Float64,
                    P_c::Float64, T_c::Float64, ρ_c::Float64, g_c::Float64,
                    P1::Float64,
                    plnt::Planet)::NamedTuple{(:L_c, :K, :ΔT, :Γ, :Ra), Float64}
    if(Ti < 0)
        P_c = plnt.P0 - plnt.a + 1e-6
        Ti = temp_adiabat(P_c, T1, P1, plnt.∇)
    end

    P_ratio = (P_c / P1)^plnt.∇

    ΔT = Ti - T_c
    ΔP = P_c - plnt.P0

    Γ = cc_slope(P_c, plnt.P0, plnt.T0, plnt.a, plnt.b)
    Γ_a = temp_distribution(P_c, T_c, plnt.∇)
    ΔΓ = Γ - Γ_a

    h = visc_height(T_c, ρ_c, g_c, plnt.A, Γ, Γ_a)
    if(h > c) h = c end

    T_m = temp_melting(P_c, plnt.P0, plnt.T0, plnt.a, plnt.b)

    η = planet_eta(plnt.η0, plnt.A, T_m, T_c)

    Ra = ρ_c * plnt.α * abs(ΔT) * g_c * abs(h)^3 / (plnt.κ * η)

    K = ρ_c * plnt.α * g_c / (plnt.Ra * η * plnt.κ)
    L_c = plnt.k * K^(1/3) * ((h / c)^(4))^(1/3) * (ΔT^4)^(1/3) * 4*π*c^2

    # this is value next to ∂c/∂t term in the derivation
    K = - 1 / (ρ_c * g_c * ΔΓ) * P_ratio

    return (L_c=L_c, K=K, ΔT=ΔT, Γ=Γ_a/ΔΓ, Ra=Ra)
end

end # module
