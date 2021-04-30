"""
    struct GadgetSphere{T} <: AbstractGadgetGeometry
        center::Vector{T}
        radius::T
    end

Defines a sphere by `center` and `radius`. 
To be used for [`read_particles_in_geometry`](@ref)
"""
struct GadgetSphere{T} <: AbstractGadgetGeometry
    center::Vector{T}
    radius::T
end


"""
    function get_geometry_center(sphere::GadgetSphere)

Returns the center of `sphere`.
"""
get_geometry_center(sphere::GadgetSphere) = sphere.center

"""
    get_geometry_box_corners(sphere::GadgetSphere)

Returns a tuple with the lower left and upper right corner of a box which contains the `sphere`.
"""
function get_geometry_box_corners(sphere::GadgetSphere)
    sphere.center .- sphere.radius, sphere.center .+ sphere.radius
end

"""
    get_geometry_mask(sphere::GadgetSphere, pos::Matrix{T}) where T

Returns the indices of all particles contained in the `sphere`.
"""
function get_geometry_mask(sphere::GadgetSphere, pos::Matrix{T}) where T
    r² = sphere.radius^2
    mask = @views @. (pos[1,:] - sphere.center[1])^2 + (pos[2,:] - sphere.center[2])^2 + (pos[3,:] - sphere.center[3])^2 ≤ r²

    return findall(mask)
end
