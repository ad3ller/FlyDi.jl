"""
    Monte-Carlo distributions
    =========================
"""

function mc_vel(num::Int64; v_s=1.95e3, delta_v=80.0, rho=2.5e-4, z=0.15, epsilon=1.0)
    """ A random sample of velocities with a mean longitudional velocity of v_s
        and a spread of delta_v.  The transverse velocity is found by assuming
        a point source and that the beam passes through a hole at z with a radius of rho
        and with x/y eccentricity of epsilon.        
    """
    # longitudional beam velocity (Gaussian)
    vz = abs.(v_s .+ randn(num) .* (delta_v / (2.0^0.5)))
    # position at the laser
    r = rho .* rand(num).^0.5
    theta = 2.0 * pi * rand(num)
    # transverse velocity (determined by skimmer)
    tof = z ./ vz
    vx = epsilon * r .* cos.(theta) ./ tof
    vy = r .* sin.(theta) ./ tof
    return Array([vx vy vz])
end

function mc_flat(num, xmin=0.0, xmax=5.0e-6)
    """ Random values, evenly sampled from xmin < x < xmax.        
    """
    x = xmin .+ rand(num) .* (xmax - xmin)
    return x
end

function mc_uniform(num::Int64; tmin=0.0, tmax=5e-6, distance=0.3, z0=0.198, reverse_vz=false,
                                v_s=1.86e3, delta_v=65.0, rho=1.5e-4, epsilon=1.0)
    """ 3D He distribution with a flat, oval distribution in xy.
    """
    t =  mc_flat(num, tmin, tmax)
    vel = mc_vel(num, v_s=v_s, delta_v=delta_v, rho=rho, z=distance, epsilon=epsilon)
    tof = distance ./ vel[:, 3]
    x = tof .* vel[:, 1]
    y = tof .* vel[:, 2]
    z = zeros(num) .+ z0
    if reverse_vz
    vel[:, 3] = - vel[:, 3]
    end
    df = DataFrame(hcat(t, x, y, z, vel))
    rename!(df, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return df
end

function mc_supersonic(num::Int64; tmin=0.0, tmax=5e-6, distance=0.2, z0=0., reverse_vz=false,
                                   v_s=2.0e3, delta_v=65.0, sigma_x=5e-4, sigma_y=1e-4, rho=1e-4)
    """ A random sample of velocities with a mean longitudional velocity of v_s
        and a spread of delta_v.  The transverse velocity is found by assuming
        a small source (rho) located at z=(z0 - distance) and that the beam has a
        2D Gaussian distribution in xy. 
    """
    # time
    t =  mc_flat(num, tmin, tmax)
    # longitudional velocity
    vz = abs.(v_s .+ randn(num) .* (delta_v / (2.0^0.5)))
    tof = distance ./ vz
    x = randn(num) .* sigma_x
    y = randn(num) .* sigma_y
    # position at the source
    r = rho .* rand(num).^0.5
    theta = 2.0 * pi * rand(num)
    x0 = r .* sin.(theta)
    y0 = r .* cos.(theta)
    # vel
    vx = (x .- x0) ./ tof
    vy = (y .- y0) ./ tof
    z = zeros(num) .+ z0
    if reverse_vz
        vz = - vz
    end
    df = DataFrame(hcat(t, x, y, z, vx, vy, vz))
    rename!(df, [:t, :x, :y, :z, :vx, :vy, :vz]);
    return df
end