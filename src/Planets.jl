# structure based off of Stixrude '21 et al.
# design still being added to it
# NOTE:
# - update T_1
# - move structure and thermal code over here
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

export Planet, Moon

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

end # module
