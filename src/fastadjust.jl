"""
    Fast adjust potential array data
    ================================

    *_r 
        indicates that function arguments include (x, y, z), i.e., real coordinates.
    
    *_g 
        indicates that function arguments include (xg, yg, zg), i.e., grid coordinates.

"""
struct FastAdjust2D
    pa::Array{Float64,3}
    el::Array{Bool,2}
    num::Int64
    nx::Int64
    ny::Int64
    dx::Float64
    dy::Float64
    x0::Float64
    y0::Float64
    xn::Float64
    yn::Float64
end;

struct FastAdjust3D
    pa::Array{Float64,4}
    el::Array{Bool,3}
    num::Int64
    nx::Int64
    ny::Int64
    nz::Int64
    dx::Float64
    dy::Float64
    dz::Float64
    x0::Float64
    y0::Float64
    z0::Float64
    xn::Float64
    yn::Float64
    zn::Float64
end;

FastAdjust = Union{FastAdjust2D, FastAdjust3D}

function h5read_pa(fil, r0::Vector{Float64}=Vector([0.0, 0.0, 0.0]))
    """ read a HDF5 file of potential array data
    """
    pa = h5read(fil, "pa")
    el = h5read(fil, "electrode")
    attrs = h5readattr(fil, "/")
    if ndims(pa) == 3
        num, ny, nx = size(pa)
        dx, dy = attrs["dx"], attrs["dy"]
        xn = r0[1] + (nx - 1) * dx
        yn = r0[2] + (ny - 1) * dy
        fa = FastAdjust2D(pa, el, num, nx, ny, dx, dy, r0[1], r0[2], xn, yn)
    elseif ndims(pa) == 4
        num, nz, ny, nx = size(pa)
        dx, dy, dz = attrs["dx"], attrs["dy"], attrs["dz"]
        xn = r0[1] + (nx - 1) * dx
        yn = r0[2] + (ny - 1) * dy
        zn = r0[3] + (nz - 1) * dz
        fa = FastAdjust3D(pa, el, num, nx, ny, nz, dx, dy, dz,
                        r0[1], r0[2], r0[3], xn, yn, zn)
    else
        ErrorException("file not compatible with FastAdjust2D or FastAdjust3D")
    end
    return fa
end

function potential(fa::FastAdjust, voltages::Vector{Float64})
    """ electric potential for a set of voltages
    """
    return squeeze(sum(fa.pa .* voltages, 1), 1)
end

function grid_r(fa::FastAdjust2D, r::Vector{Float64})
    """ convert real coordinates into grid coordinates
    """
    xg = 1.0 + (r[1] - fa.x0) / fa.dx
    yg = 1.0 + (r[2] - fa.y0) / fa.dy
    zg = NaN
    return xg, yg, zg
end

function grid_r(fa::FastAdjust3D, r::Vector{Float64})
    """ convert real coordinates into grid coordinates
    """
    xg = 1.0 + (r[1] - fa.x0) / fa.dx
    yg = 1.0 + (r[2] - fa.y0) / fa.dy
    zg = 1.0 + (r[3] - fa.z0) / fa.dz
    return xg, yg, zg
end

function electrode_g(fa::FastAdjust2D, xg::Float64, yg::Float64, zg::Float64)
    """ are any of the grid point adjacent to (xg, yg, zg) an electrode?
    """
    xn = Int64(floor(xg))
    yn = Int64(floor(yg))
    if (1 <= xn < fa.nx) && (1 <= yn < fa.ny)
        return any(fa.el[yn : yn + 1, xn : xn + 1])
    else
        # everywhere outside the pa is an electrode
        return true
    end
end

function electrode_g(fa::FastAdjust3D, xg::Float64, yg::Float64, zg::Float64)
    """ are any of the grid point adjacent to (xg, yg, zg) an electrode?
    """
    xn = Int64(floor(xg))
    yn = Int64(floor(yg))
    zn = Int64(floor(zg))
    if (1 <= xn < fa.nx) && (1 <= yn < fa.ny) && (1 <= zn < fa.nz)
        return any(fa.el[zn : zn + 1, yn : yn + 1, xn : xn + 1])
    else
        # everywhere outside the pa is an electrode
        return true
    end
end

function potential_g(fa::FastAdjust2D, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ electric potential at grid coord
    """
    try
        # integer coords (floor rounding)
        xn = Int64(floor(xg))
        yn = Int64(floor(yg))
        # get potential
        wx = xg - xn
        wy = yg - yn
        pa = fa.pa[:, yn: yn + 1, xn : xn + 1]
        phi = squeeze(sum(pa .* voltages, 1), 1)
        result = (1.0 - wx) * (1.0 - wy) * phi[1, 1] + 
                 wx * (1.0 - wy) * phi[1, 2] + 
                 (1.0 - wx) * wy * phi[2, 1] +
                 wx * wy * phi[2, 2]
        return result   
    catch
        return NaN
    end
end

function potential_g(fa::FastAdjust3D, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ electric potential at grid coord
    """
    try
        # integer coords (floor rounding)
        xn = Int64(floor(xg))
        yn = Int64(floor(yg))
        zn = Int64(floor(zg))
        pa = fa.pa[:, zn: zn + 1, yn: yn + 1, xn : xn + 1]
        phi = squeeze(sum(pa .* voltages, 1), 1)
        ## interpolate along x
        wx = (xg - xn)
        c00 = phi[1, 1, 1] * (1 - wx) .+ phi[1, 1, 2] * wx
        c01 = phi[2, 1, 1] * (1 - wx) .+ phi[2, 1, 2] * wx
        c10 = phi[1, 2, 1] * (1 - wx) .+ phi[1, 2, 2] * wx
        c11 = phi[2, 2, 1] * (1 - wx) .+ phi[2, 2, 2] * wx
        ## interpolate along y
        wy = (yg - yn)
        c0 = c00 * (1 - wy) .+ c10 * wy
        c1 = c01 * (1 - wy) .+ c11 * wy
        ## interpolate along z
        wz = (zg - zn)
        return c0 * (1 - wz) .+ c1 * wz
    catch
        return NaN
    end
end

function field_g(fa::FastAdjust2D, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ electric field at grid coord
    """
    ex = (potential_g(fa, voltages, xg + 0.5, yg, zg) - 
          potential_g(fa, voltages, xg - 0.5, yg, zg)) / fa.dx
    ey = (potential_g(fa, voltages, xg, yg + 0.5, zg) - 
          potential_g(fa, voltages, xg, yg - 0.5, zg)) / fa.dy
    ez = 0.0
    return ex, ey, ez
end

function field_g(fa::FastAdjust3D, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ electric field at grid coord
    """
    ex = (potential_g(fa, voltages, xg + 0.5, yg, zg) - 
          potential_g(fa, voltages, xg - 0.5, yg, zg)) / fa.dx
    ey = (potential_g(fa, voltages, xg, yg + 0.5, zg) - 
          potential_g(fa, voltages, xg, yg - 0.5, zg)) / fa.dy
    ez = (potential_g(fa, voltages, xg, yg, zg + 0.5) - 
          potential_g(fa, voltages, xg, yg, zg - 0.5)) / fa.dz
    return ex, ey, ez
end

function amp_field_g(fa::FastAdjust, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ amplitude of the field at grid coord
    """
    ex, ey, ez = field_g(fa, voltages, xg, yg, zg)
    return sqrt(ex^2.0 + ey^2.0 + + ez^2.0)
end

function grad_field_g(fa::FastAdjust2D, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ gradient in the amplitude of the electric field at grid coord
    """
    dfx = (amp_field_g(fa, voltages, xg + 0.5, yg, zg) - 
           amp_field_g(fa, voltages, xg - 0.5, yg, zg)) / fa.dx
    dfy = (amp_field_g(fa, voltages, xg, yg + 0.5, zg) - 
           amp_field_g(fa, voltages, xg, yg - 0.5, zg)) / fa.dy
    dfz = 0.0
    return dfx, dfy, dfz
end

function grad_field_g(fa::FastAdjust3D, voltages::Vector{Float64}, xg::Float64, yg::Float64, zg::Float64)
    """ gradient in the amplitude of the electric field at grid coord
    """
    dfx = (amp_field_g(fa, voltages, xg + 0.5, yg, zg) - 
           amp_field_g(fa, voltages, xg - 0.5, yg, zg)) / fa.dx
    dfy = (amp_field_g(fa, voltages, xg, yg + 0.5, zg) - 
           amp_field_g(fa, voltages, xg, yg - 0.5, zg)) / fa.dy
    dfz = (amp_field_g(fa, voltages, xg, yg, zg + 0.5) - 
           amp_field_g(fa, voltages, xg, yg, zg - 0.5)) / fa.dz
    return dfx, dfy, dfz
end
