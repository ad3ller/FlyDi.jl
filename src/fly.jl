"""
    Fly dipoles
    ===========

    flydi
        trajectories for a range of starting conditions
    
    flydi_last 
        final position for a range of starting conditions

    Notes
    -----
    The function results are in the same format as the starting input, i.e., Array -> Array, DataFrame -> DataFrame.

"""

function flydi(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
               di::DataFrame,
               dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
               max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ M-C trajectory (DataFrames)
    """
    # fly dipoles
    num, cols = size(di)
    result = Array{DataFrame}(num)
    for i in 1:num
        t0, r0, v0 = di[i, :t], Array(di[i, [:x, :y, :z]])[1, :], Array(di[i, [:vx, :vy, :vz]])[1, :]
        result[i] = traj(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                         max_field=max_field, max_t=max_t, max_iterations=max_iterations, df=true)
    end
    return result
end

function flydi(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
               t_initial::Vector{Float64}, r_initial::Array{Float64, 2}, v_initial::Array{Float64, 2},
               dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
               max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ M-C trajectory (Array)
    """
    # fly dipoles
    num, dims = size(r_initial)
    result = Array{Array{Float64, 2}}(num)
    for i in 1:num
        t0, r0, v0 = t_initial[i], r_initial[i, :], v_initial[i, :]
        result[i] = traj(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                         max_field=max_field, max_t=max_t, max_iterations=max_iterations, df=false)
    end
    return result
end

function flydi_last(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                    di::DataFrame,
                    dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                    max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6), parallel::Bool=false)
    """ M-C final position (DataFrame)
    """
    # fly dipoles
    t_initial = Array(di[:t])
    r_initial = Array(di[[:x, :y, :z]])
    v_initial = Array(di[[:vx, :vy, :vz]])
    result = flydi_last(fa, voltages, t_initial, r_initial, v_initial, dipole, mass, dt, method,
                        max_field=max_field, max_t=max_t, max_iterations=max_iterations, parallel=parallel)           
    result = DataFrame(result)
    names!(result, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return result
end

function flydi_last(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                    t_initial::Vector{Float64}, r_initial::Array{Float64, 2}, v_initial::Array{Float64, 2},
                    dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                    max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6), parallel::Bool=false)
    """ M-C final position (Array)
    """
    if parallel
        result = flydi_last_parallel(fa, voltages, t_initial, r_initial, v_initial, dipole, mass, dt, method,
                              max_field=max_field, max_t=max_t, max_iterations=max_iterations)
    else
        result = flydi_last_serial(fa, voltages, t_initial, r_initial, v_initial, dipole, mass, dt, method,
                              max_field=max_field, max_t=max_t, max_iterations=max_iterations)
    end
    return result
end

function flydi_last_serial(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                    t_initial::Vector{Float64}, r_initial::Array{Float64, 2}, v_initial::Array{Float64, 2},
                    dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                    max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ M-C final position (Array)
    """
    # fly dipoles
    num, dims = size(r_initial)
    result = Array{Float64}(num, 7)
    for i in 1:num
        t0, r0, v0 = t_initial[i], r_initial[i, :], v_initial[i, :]
        result[i, :] = traj_last(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                                 max_field=max_field, max_t=max_t, max_iterations=max_iterations)
    end
    return result
end

function flydi_last_parallel(fa::FastAdjust, voltages::Union{Function, Vector{Float64}},
                      t_initial::Vector{Float64}, r_initial::Array{Float64, 2}, v_initial::Array{Float64, 2},
                      dipole::Union{Float64, Function}, mass::Float64, dt::Float64, method::Function;
                      max_field::Float64=NaN, max_t::Float64=1.0, max_iterations::Int64=Int64(1e6))
    """ M-C final position (parallel)
    """
    num, dims = size(r_initial)
    result = SharedArray{Float64}(num, 7)
    @sync @parallel for i in 1:num
        t0, r0, v0 = t_initial[i], r_initial[i, :], v_initial[i, :]
        result[i, :] = traj_last(fa, voltages, t0, r0, v0, dipole, mass, dt, method,
                                 max_field=max_field, max_t=max_t, max_iterations=max_iterations);
    end
    return result
end