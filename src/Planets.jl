# structure based off of Stixrude '21 et al.
# design still being added to it
module Planets

export Planet

struct Planet

    # structure parameters
    layers :: Real        # number of layers in structure file
    R      :: Real        # radius of planet
    M      :: Real        # mass of planet

    # orbital parameters
    ω :: Real           # rotational frequency

    # thermal parameters
    C_p  :: Real        # specific heat
    α    :: Real        # thermal expansivity
    k    :: Real        # thermal conductivity
    T0   :: Real        # phase transition reference T
    P0   :: Real        # phase transition reference P
    a    :: Real        # Simon pressure
    b    :: Real        # Simon gradient
    B    :: Real        # effective temperature exponent
    ∇    :: Real        # adiabtatic gradient
    T_eq :: Real        # radiative equilibrium temperature

    # rheaology parameters
    η0         :: Real    # reference viscosity
    A          :: Real    # viscosity exponent
    Ra         :: Real    # critical Rayleigh number
    γ          :: Real    # Andrade exponent
    μ_f        :: Real    # standard linear solid parameter
    rhea_model :: Real    # modulus model to use
                            # 1 - Maxwell, 2 - SLS, 3 - Andrade

end

end # module
