# structure based off of Stixrude '21 et al.
# design still being added to it
module Planets

export Planet, Moon

struct Planet

    name :: String     # name of planet

    # structure parameters
    layers :: Real        # number of layers in structure file
    R      :: Real        # radius of planet
    M      :: Real        # mass of planet
    GM     :: Real        # gravitational mass of planet

    # orbital parameters
    ω :: Real           # rotational frequency

    # thermal parameters
    C_p  :: Real        # specific heat
    α    :: Real        # thermal expansivity
    k    :: Real        # thermal conductivity
    κ    :: Real        # thermal diffusivity
    T0   :: Real        # phase transition reference T
    P0   :: Real        # phase transition reference P
    a    :: Real        # Simon pressure
    b    :: Real        # Simon gradient
    B    :: Real        # effective temperature exponent
    ∇    :: Real        # adiabtatic gradient
    T_eq :: Real        # radiative equilibrium temperature
    T_ef :: Real        # Today's effective temperature
    T1   :: Real        # Temp at 1 bar

    # rheaology parameters
    η0         :: Real    # reference viscosity
    A          :: Real    # viscosity exponent
    Ra         :: Real    # critical Rayleigh number
    μ_f        :: Tuple   # rhealogoy parameters
                            # 1 - SLS mu factor, 2,3 - Andrade alpha, beta
    rhea_model :: Real    # modulus model to use
                            # 1 - Maxwell, 2 - SLS, 3 - Andrade

end

struct Moon

    name :: String    # name of moon

    # structure parameters
    m  :: Real        # mass of moon
    gm :: Real        # gravitational mass of moon

    # orbital parameters
    a :: Real         # orbital radius
    i :: Real         # inclination

    # flags
    retro :: Bool     # is planet in retrograde motion

end

end # module
