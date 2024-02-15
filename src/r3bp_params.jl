"""
Function define CR3BP system parameters
"""

mutable struct CR3BP_param
    mu::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    mstar::Float64
    m2_soi::Float64
end


"""
    get_cr3bp_param(m1_naifID::Int, m2_naifID::Int)

Obtain CR3BP parameters mu, Lstar, Tstar, soi of m2

Args:
    m1_naifID (Int): mass of first primary
    m2_naifID (Int): mass of second primary

Returns:
    (obj): object with fields:
        mu (float): mass-parameter
        lstar (float): non-dimensional distance
        tstar (float): non-dimensional time
        m2_soi (float): sphere of influence of second mass, in km
"""
function get_cr3bp_param(m1_naifID::Int, m2_naifID::Int)
    return get_cr3bp_param(string(m1_naifID), string(m2_naifID))
end



"""
    get_cr3bp_param(m1_naifID::String, m2_naifID::String)

Obtain CR3BP parameters mu, lstar, tstar, soi of m2

# Arguments
    m1_naifID (str): mass of first primary
    m2_naifID (str): mass of second primary

# Returns
    (obj): object with fields:
        mu (float): mass-parameter
        lstar (float): non-dimensional distance
        tstar (float): non-dimensional time
        m2_soi (float): sphere of influence of second mass, in km
"""
function get_cr3bp_param(m1_naifID::String, m2_naifID::String)
    # list of semi-major axis
    if m1_naifID =="10"
        a2 = get_semiMajorAxes(m2_naifID)[1]
    elseif m2_naifID =="10"
        a2 = get_semiMajorAxes(m1_naifID)[1]
    else
        semiMajorAxes = get_semiMajorAxes(m1_naifID, m2_naifID)
        a1, a2 = semiMajorAxes[1], semiMajorAxes[2]
    end
    # list of gm
    gmlst = get_gm(m1_naifID, m2_naifID)
    m1_gm, m2_gm = gmlst[1], gmlst[2]
    if m1_gm < m2_gm
        error("Excepting m1 > m2!")
    end
    # create list of parameters
    mu     = m2_gm / (m1_gm + m2_gm)
    lstar  = a2
    tstar  = sqrt( ( a2 )^3 / ( m1_gm + m2_gm ) )
    vstar = lstar/tstar
    mstar  = m1_gm + m2_gm
    m2_soi = a2 * (m2_gm/m1_gm)^(2/5)
    return CR3BP_param(mu, lstar, tstar, vstar, mstar, m2_soi)
end


mutable struct BCR4BP_param
    mu::Float64
    lstar::Float64
    tstar::Float64
    vstar::Float64
    mstar::Float64
    m2_soi::Float64
    μ_3::Float64
    a::Float64
    ω_s::Float64
    tsyn::Float64
end


"""
    get_bcr4bp_param(m1_naifID::Int, m2_naifID::Int)

Obtain BCR4BP parameters mu, Lstar, Tstar, soi of m2, from m1-m2 planet-moon-sun system

Args:
    m1_naifID (Int): mass of first primary
    m2_naifID (Int): mass of second primary

Returns:
    (obj): object with fields:
        mu (float): mass-parameter
        lstar (float): non-dimensional distance
        tstar (float): non-dimensional time
        m2_soi (float): sphere of influence of second mass, in km
"""
function get_bcr4bp_param(m1_naifID::Int, m2_naifID::Int)
    return get_bcr4bp_param(string(m1_naifID), string(m2_naifID))
end


"""
    get_bcr4BP_param(m1_naifID::String, m2_naifID::String)

Obtain BCR4BP parameters mu, lstar, tstar, soi of m2, from m1-m2 planet-moon-sun system

Args:
    m1_naifID (str): mass of first primary
    m2_naifID (str): mass of second primary

Returns:
    (obj): object with fields:
        mu (float): mass-parameter
        lstar (float): non-dimensional distance
        tstar (float): non-dimensional time
        m2_soi (float): sphere of influence of second mass, in km
"""
function get_bcr4bp_param(m1_naifID::String, m2_naifID::String)
    CR3BP_param = get_cr3bp_param(m1_naifID, m2_naifID)

    # get information of Sun
    gm_sun     = get_gm("10")[1]
    a_m1_sun   = get_semiMajorAxes(m1_naifID)[1]   # semi-major axis of first (larger) body
    println("a_m1_sun: $a_m1_sun, gm_sun: $gm_sun")
    period_sun = 2π*sqrt(a_m1_sun^3/gm_sun)   # [sec]

    # compute rotation rate of sun about m1-m2 barycenter
    tsyn = get_synodic_period(period_sun, 2π*CR3BP_param.tstar)  # [sec]
    ω_s  = -2π / (tsyn / CR3BP_param.tstar)                    # [rad/canonical time]
    println("ω_s: $ω_s")
    # compute scaled gm and length of sun
    a_sun = a_m1_sun / CR3BP_param.lstar
    μ_3   = gm_sun   / CR3BP_param.mstar

    # return structure
    return BCR4BP_param(CR3BP_param.mu, CR3BP_param.lstar, CR3BP_param.tstar, CR3BP_param.vstar,
    CR3BP_param.mstar, CR3BP_param.m2_soi, μ_3, a_sun, ω_s, tsyn)
end



"""
    get_synodic_period(p1::Real, p2::Real)

Compute synodic period between two systems with periods p1 and p2

# Arguments
    p1 (float): period of first system
    p2 (float): period of second system

# Returns
    (float): synodic period
"""
function get_synodic_period(p1::Real, p2::Real)
    return 1/abs( 1/p1 - 1/p2 )
end


"""
    get_resonance_period(M::Real, N::Real)

Compute M-N resonance period
"""
function get_resonance_period(M::Real, N::Real)
    return 2*N*π/M
end


"""
    get_jacobi_constant(μ::Float64, sv)

Compute Jacobi constant

# Arguments
    - `μ::Float64`: CR3BP parameter
    - `sv::Vector`: state-vector, planar or spatial
"""
function jacobi_constant(μ::Float64, sv)
    if length(sv)==4
        x, y, vx, vy = sv
        z, vz = 0.0, 0.0;
    elseif length(sv)==6
        x, y, z, vx, vy, vz = sv
    end
    r1 = sqrt( (x+μ)^2 + y^2 + z^2 );
    r2 = sqrt( (x-1+μ)^2 + y^2 + z^2 );
    μ1 = 1.0-μ;
    μ2 = μ;
    # compute augmented potential
    ubar = 0.5*(x^2 + y^2) + μ1/r1 + μ2/r2;
    jc   = 2*ubar - (vx^2 + vy^2 + vz^2);
    return jc
end


"""
    lagrange_points(μ::Float64)

Function computes lagrange points from CR3BP paraneter μ"""
function lagrange_points(μ::Float64)
    # place-holder
    l = 1 - μ;
    # collinear points
    fl1(x) = x^5 + 2*(μ-l)*x^4 + (l^2-4*l*μ+μ^2)*x^3 + (2*μ*l*(l-μ)+μ-l)*x^2 + (μ^2*l^2+2*(l^2+μ^2))*x + (μ^3-l^3);
    xl1 = Roots.find_zero(fl1, (0,2), Roots.Bisection());

    fl2(x) = x^5 + 2*(μ-l)*x^4 + (l^2-4*l*μ+μ^2)*x^3 + (2*μ*l*(l-μ)-μ-l)*x^2 + (μ^2*l^2+2*(l^2-μ^2))*x - (μ^3+l^3);
    xl2 = Roots.find_zero(fl2, (0,2), Roots.Bisection());

    fl3(x) = x^5 + 2*(μ-l)*x^4 + (l^2-4*l*μ+μ^2)*x^3 + (2*μ*l*(l-μ)+μ+l)*x^2 + (μ^2*l^2+2*(μ^2-l^2))*x + (μ^3+l^3);
    xl3 = Roots.find_zero(fl3, (-2,0), Roots.Bisection());

    LP = zeros(5, 6);
    LP[1,1] = xl1;
    LP[2,1] = xl2;
    LP[3,1] = xl3;
    # equilateral points
    LP[4,1] = cos(pi/3)-μ;
    LP[4,2] = sin(pi/3);
    LP[5,1] = cos(pi/3)-μ;
    LP[5,2] = -sin(pi/3);
    return LP
end
