"""
    Fly dipoles
    ===========

    flydi
        trajectories for a range of starting conditions
    
    flydi_last 
        final position for a range of starting conditions

"""

function flydi(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
               di::DataFrame,
               dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
               max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ Dipole trajectory
    """
    # fly dipoles
    num, cols = size(di)
    result = Array{DataFrame, 1}(undef, num)
    for i in 1:num
        t0, r0, v0 = di[i, :t], Vector(di[i, [:x, :y, :z]]), Vector(di[i, [:vx, :vy, :vz]])
        result[i] = traj(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                         max_field=max_field, max_t=max_t, max_iterations=max_iterations, df=true)
    end
    return result
end

function flydi_last(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                     di::DataFrame,
                     dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                     max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ Dipole trajectory (final position)
    """
    # fly dipoles
    t_initial = Array(di[:t])
    r_initial = convert(Matrix, di[[:x, :y, :z]])
    v_initial = convert(Matrix, di[[:vx, :vy, :vz]])
    num, dims = size(r_initial)
    result = Array{Float64, 2}(undef, num, 7)
    for i in 1:num
        t0, r0, v0 = t_initial[i], r_initial[i, :], v_initial[i, :]
        result[i, :] = traj_last(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                                 max_field=max_field, max_t=max_t, max_iterations=max_iterations)
    end
    result = DataFrame(result)
    rename!(result, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return result
end
