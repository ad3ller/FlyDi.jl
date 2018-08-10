module FlyDi

export c, q, u, m_e, h, a0, Ry, kB,
    FastAdjust,
    h5read_pa, potential,
    grid_r, electrode_g, pa_g, potential_g, field_g, amp_field_g, grad_field_g,
    kinetic_energy, potential_energy,
    euler, leapfrog, rk4,
    traj, traj_last,
    flydi, flydi_last,
    mc_supersonic, mc_uniform

include("external.jl")
include("constants.jl")
include("fastadjust.jl")
include("traj.jl")
include("fly.jl")
include("montecarlo.jl")

end