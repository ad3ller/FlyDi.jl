"""
    Trajectory of a dipole in an electric field
    ===========================================
"""

function voltages_t(voltages::Vector{Float64}, t::Float64)
    """ electrode voltages at time t
    """
    return voltages
end

function voltages_t(voltages::Function, t::Float64)
    """ electrode voltages at time t
    """
    return voltages(t)
end

function kinetic_energy(v::Vector{Float64}, mass::Float64)
    """ kinetic energy of a moving particle
    """
    return 0.5 * mass * (v[1]^2.0 + v[2]^2.0 + v[3]^2.0)
end

function potential_energy(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, xg::Float64, yg::Float64, zg::Float64, dipole::Float64; get_famp::Bool=false)
    """ potential energy of a static dipole
    """
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)
    if get_famp
        return dipole * famp, famp
    else
        return dipole * famp
    end
end

function potential_energy(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, xg::Float64, yg::Float64, zg::Float64, dipole::Function; get_famp::Bool=false)
    """ potential energy of a field-dependent dipole
    """
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)
    dpm = dipole(famp)
    if get_famp
        return dpm * famp, famp
    else
        return dpm * famp
    end
end

function force(fa::FastAdjust3D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole_f::Function)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (3D)
    """
    xg, yg, zg = grid_r(fa, r)
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)
    dipole = dipole_f(famp)
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy
    fz = -(potential_energy(fa, voltages, t, xg, yg, zg + 0.5, dipole) - 
           potential_energy(fa, voltages, t, xg, yg, zg - 0.5, dipole)) / fa.dz
    return Vector([fx, fy, fz])
end

function force(fa::FastAdjust2D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole_f::Function)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (2D)
    """
    xg, yg, zg = grid_r(fa, r)
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)
    dipole = dipole_f(famp)
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy
    fz = 0.0
    return Vector([fx, fy, fz])
end

function force(fa::FastAdjust3D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole::Float64)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (3D)
    """
    xg, yg, zg = grid_r(fa, r)
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy
    fz = -(potential_energy(fa, voltages, t, xg, yg, zg + 0.5, dipole) - 
           potential_energy(fa, voltages, t, xg, yg, zg - 0.5, dipole)) / fa.dz
    return Vector([fx, fy, fz])
end

function force(fa::FastAdjust2D, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole::Float64)
    """ gradient of the potential energy of a dipole in an inhomegeous electric field (2D)
    """
    xg, yg, zg = grid_r(fa, r)
    fx = -(potential_energy(fa, voltages, t, xg + 0.5, yg, zg, dipole) - 
           potential_energy(fa, voltages, t, xg - 0.5, yg, zg, dipole)) / fa.dx
    fy = -(potential_energy(fa, voltages, t, xg, yg + 0.5, zg, dipole) - 
           potential_energy(fa, voltages, t, xg, yg - 0.5, zg, dipole)) / fa.dy
    fz = 0.0
    return Vector([fx, fy, fz])
end

function acceleration(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64)
    """ accelleration of a dipole in an inhomegeous electric field
    """
    return  force(fa, voltages, t, r, dipole) ./ mass
end

function euler(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, v::Vector{Float64}, a::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64, dt::Float64)
    """ Euler method of integration
    """
    a = acceleration(fa, voltages, t, r, dipole, mass)
    rn = r .+ v .* dt
    vn = v .+ a .* dt
    return rn, vn, a
end

function leapfrog(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, v::Vector{Float64}, a::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64, dt::Float64)
    """ leapfrog method of integration
    """
    if any(isnan.(a))
        a = acceleration(fa, voltages, t, r, dipole, mass)
    end
    rn = r .+ v .* dt + a .* (0.5 * dt^2.0)
    an = acceleration(fa, voltages, t + dt, rn, dipole, mass)
    vn = v .+ (a .+ an) .* (0.5 * dt)
    return rn, vn, an
end

function rk4(fa::FastAdjust, voltages::Union{Function, Vector{Float64}}, t::Float64, r::Vector{Float64}, v::Vector{Float64}, a::Vector{Float64}, dipole::Union{Float64, Function}, mass::Float64, dt::Float64)
    """ fourth order Runge-Kutta method of integration
    """
    # initial position
    a = acceleration(fa, voltages, t, r, dipole, mass)
    # half-step
    r1 = r .+ (0.5 * dt) .* v
    v1 = v .+ (0.5 * dt) .* a
    a1 = acceleration(fa, voltages, t + 0.5 * dt, r1, dipole, mass)
    # half-step again
    r2 = r .+ (0.5 * dt) .* v1
    v2 = v .+ (0.5 * dt) .* a1
    a2 = acceleration(fa, voltages, t + 0.5 * dt, r2, dipole, mass)
    # full step
    r3 = r .+ dt .* v2
    v3 = v .+ dt .* a2
    a3 = acceleration(fa, voltages, t + dt, r3, dipole, mass)
    # next position
    r4 = r .+ (dt / 6.0) .* (v .+ 2.0 .* v1 .+ 2.0 .* v2 .+ v3)
    v4 = v .+ (dt / 6.0) .* (a .+ 2.0 .* a1 .+ 2.0 .* a2 .+ a3)
    return r4, v4, a3
end

function traj(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
              t0::Float64, r0::Vector{Float64}, v0::Vector{Float64},
              dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
              max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6), df::Bool=false)
    """ full trajectory
    """
    # initialise
    i = 1
    t, r, v = t0, r0, v0
    a = Vector([NaN, NaN, NaN])
    xg, yg, zg = grid_r(fa, r)
    KE = kinetic_energy(v, mass)
    PE = potential_energy(fa, voltages, t, xg, yg, zg, dipole)
    famp = amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg)
    result = [t r[1] r[2] r[3] KE PE famp;]
    # step-by-step trajectory
    while i < max_iterations
        try
            r, v, a = method(fa, voltages, t, r, v, a, dipole, mass, dt)
            t += dt
        catch
            break
        end
        if any(isnan, r) | any(isnan, v)
            break
        end
        # grid position
        xg, yg, zg = grid_r(fa, r)
        # energy
        KE = kinetic_energy(v, mass)
        PE, famp = potential_energy(fa, voltages, t, xg, yg, zg, dipole, get_famp=true)
        # checks
        if electrode_g(fa, xg, yg, zg)
            # hit an electrode or left the pa
            break
        elseif ~isnan(max_field) && famp > max_field
            # ionisation field
            break
        elseif ~isnan(max_t) && t > max_t
            # time's up!
            break
        end
        # record
        result = vcat(result, [t r[1] r[2] r[3] KE PE famp])
        # next step
        i += 1
    end
    # output
    if df
        result = DataFrame(result)
        names!(result, [:t, :x, :y, :z, :KE, :PE, :famp]);
    end
    return result
end

function traj_last(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                   t0::Float64, r0::Vector{Float64}, v0::Vector{Float64},
                   dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                   max_field::Float64=NaN, max_t::Float64=NaN, max_iterations::Int64=Int64(1e6))
    """ last point of a trajectory
    """
    # initialise
    i = 1
    t, r, v = t0, r0, v0
    a = Vector([NaN, NaN, NaN])
    result = [t r[1] r[2] r[3] v[1] v[2] v[3];]
    # step-by-step trajectory
    while i < max_iterations
        try
            r, v, a = method(fa, voltages, t, r, v, a, dipole, mass, dt)
            t += dt
        catch
            break
        end
        if any(isnan, r) | any(isnan, v)
            break
        end
        # grid position
        xg, yg, zg = grid_r(fa, r)
        # checks
        if electrode_g(fa, xg, yg, zg)
            # hit an electrode or left the pa
            break
        elseif ~isnan(max_field) && amp_field_g(fa, voltages_t(voltages, t), xg, yg, zg) > max_field
            # ionisation field
            break
        elseif ~isnan(max_t) && t > max_t
            # time's up!
            break
        end
        # record
        result = [t r[1] r[2] r[3] v[1] v[2] v[3];]
        # next step
        i += 1
    end
    # output
    return result
end
