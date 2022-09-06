"""
    struct GadgetShell{T} <: AbstractGadgetGeometry
        center::Vector{T}
        radius::T
        width::T
    end

Defines a shell by `center`, `radius` and `width`.
Radius is defined to the middle of the shell thickness. 
To be used for [`read_particles_in_geometry`](@ref)
"""
struct GadgetShell{T} <: AbstractGadgetGeometry
    center::Vector{T}
    radius::T
    width::T
end


"""
    function get_geometry_center(shell::GadgetShell)

Returns the center of `shell`.
"""
get_geometry_center(shell::GadgetShell) = shell.center

"""
    get_geometry_box_corners(shell::GadgetShell)

Returns a tuple with the lower left and upper right corner of a box which contains the `shell`.
"""
function get_geometry_box_corners(shell::GadgetShell)
    dx = shell.radius + 0.5shell.width
    shell.center .- dx, shell.center .+ dx
end

"""
    get_geometry_mask(shell::GadgetShell, pos::Matrix{T}) where T

Returns the indices of all particles contained in the `shell`.
"""
function get_geometry_mask(shell::GadgetShell, pos::Matrix{T}) where T
    r1²  = (shell.radius + 0.5shell.width)^2
    r2²  = (shell.radius - 0.5shell.width)^2

    d²   = @views @. (pos[1,:] - shell.center[1])^2 + (pos[2,:] - shell.center[2])^2 + (pos[3,:] - shell.center[3])^2
    mask = @views @. ( d² ≤ r1² ) & ( d² ≥ r2² )

    return findall(mask)
end
