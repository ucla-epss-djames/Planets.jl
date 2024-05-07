export planet_eta
export temp_effective, temp_adiabat, temp_melting
export cc_slope, temp_distribution, visc_height, thermal_diff, layer_density

"""
    planet_eta(η0::Float64, A::Float64, T_m::Float64, T::Float64)::Float64

Calculates viscosity. Refer to Stixrude et al. 2021 (eq 15).

# Arguments
- `η0::Float64`  - reference viscosity
- `A::Float64`   - viscosity exponent
- `T_m::Float64` - melting temperature
- `T::Float64`   - temperature
"""
function planet_eta(η0::Float64, A::Float64, T_m::Float64, T::Float64)::Float64
    η0 * exp(A * (T_m/T - 1))
end

"""
    temp_effective(T1::Union{Int64, Float64}, B::Float64)::Float64

Calcultes the effective temperature where it is the temperature the planet would
have in absence of solar luminosity. Refer to Stixrude et al. 2021 (eq 13) or
Guillot et al. 1995.

# Arguments
- `T1::Union{Int64, Float64}` - temperature
- `B::Float64`  - constant of normalization
"""
function temp_effective(T1::Union{Int64, Float64}, B::Float64)::Float64
    (T1 / B)^(1/1.244)
end

"""
    temp_adiabat(P::Float64, Tx::Float64, Px::Float64, ∇::Float64)::Float64

Calculates the adiabatic temperature. Refer to Stixrude et al. 2021 (eq 12).

# Arguments
- `P::Float64`  - pressure
- `Tx::Float64` - temperature
- `Px::Float64` - reference pressure
- `∇::Floa64`  - adiabatic gradient
"""
function temp_adiabat(P::Float64, Tx::Float64, Px::Float64, ∇::Float64)::Float64
    Tx * (P / Px)^∇
end

"""
    temp_melting(P::Float64, P0::Union{Int64, Float64},
                 T0::Union{Int64, Float64}, a::Float64, b::Float64)::Float64

Calculates the melting temperature. Refer to Stixrude et al. 2021 (eq 11).

# Arguments
- `P::Float64`   - pressure
- `P0::Union{Int64, Float64}` - phase transition reference pressure
- `T0::Union{Int64, Float64}` - phase transition reference temperature
- `a::Float64`   - Simon pressure
- `b::Float64`   - Simon exponent
"""
function temp_melting(P::Float64, P0::Union{Int64, Float64},
                      T0::Union{Int64, Float64}, a::Float64, b::Float64)::Float64
    T0 * (1 + (P - P0) / a)^b
end

"""
    cc_slope(P::Float64, P0::Union{Int64, Float64},
             T0::Union{Int64, Float64}, a::Float64, b::Float64)::Float64

Calculates the Clausius-Clapeyron slope of the phase boundary, Γ, where
Γ = d(T_m)/dP.

# Arguments
- `P::Float64`   - pressure
- `P0::Union{Int64, Float64}`  - phase transition reference pressure
- `T0::Union{Int64, Float64}`  - phase transition reference temperature
- `a::Float64`   - Simon pressure
- `b::Float64`   - Simon exponent
"""
function cc_slope(P::Float64, P0::Union{Int64, Float64},
                  T0::Union{Int64, Float64}, a::Float64, b::Float64)::Float64
    T0 * b / a * (1 + (P - P0) / a)^(b - 1)
end

"""
    temp_distribution(P_c::Float64, T_c::Float64, ∇::Float64)::Float64

Calculates the slope, Γ_a, of the planetary temperature distribution at `T_c` and
`P_c`.

# Arguments
- `P_c::Float64` - pressure at top of the core
- `T_c::Float64` - temperature at top of the core
- `∇::Float64`   - adiabatic gradient
"""
function temp_distribution(P_c::Float64, T_c::Float64, ∇::Float64)::Float64
    (T_c / P_c) * ∇
end

"""
    visc_height(T_c::Float64, ρ_c::Float64, g_c::Float64, A::Int64,
                Γ::Float64, Γ_a::Float64)::Float64

Calculates the viscous scale height, h.

# Arguments
- `T_c::Float64` - temperature at top of the core
- `ρ_c::Float64` - density at top of the core
- `g_c::Float64` - gravity at top of the core
- `A::Int64`     - viscosity exponent
- `Γ::Float64`   - Clausius-Clapeyron slope of phase boundary
- `Γ_a::Float64` - slope of planetary temp distribution
"""
function visc_height(T_c::Float64, ρ_c::Float64, g_c::Float64, A::Int64,
                     Γ::Float64, Γ_a::Float64)::Float64
    T_c / (ρ_c * g_c * A * (Γ - Γ_a))
end

"""
    thermal_diff(k::Float64, ρ::Float64, C_p::Union{Int64, Float64})::Float64

Calculates the thermal diffusivity.

# Arguments
- `k::Float64`   - thermal conductivity
- `ρ::Float64`   - density
- `C_p::Union{Int64, Float64}` - specific heat
"""
function thermal_diff(k::Float64, ρ::Float64, C_p::Union{Int64, Float64})::Float64
    k / (ρ * C_p)
end

"""
    layer_density(r::Float64, P::Float64, Px::Float64, ∇::Float64,
                  ρ::Float64)::Float64

Function integrated for luminosity. Refer to Stixrude et al. 2021 (eq 8-9) or
Fortney & Nettelmann 2010.

# Arguments
- `r::Float64`  - radius
- `P::Float64`  - pressure
- `Px::Float64` - reference pressure
- `∇::Float64`  - adiabatic gradient
- `ρ::Float64`  - density
"""
function layer_density(r::Float64, P::Float64, Px::Float64, ∇::Float64,
                       ρ::Float64)::Float64
    (P/Px)^∇ * ρ * r^2
end
