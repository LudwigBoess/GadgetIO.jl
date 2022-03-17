"""
    struct GadgetCylinder{T} <: AbstractGadgetGeometry
        pos_start::Vector{T}
        pos_end::Vector{T}
        radius::T
    end

Defines a cylinder by two end points `pos_start` and `pos_start` as well as a `radius`. 
To be used for [`read_particles_in_geometry`](@ref)
"""
struct GadgetCylinder{T} <: AbstractGadgetGeometry
    pos_start::Vector{T}
    pos_end::Vector{T}
    radius::T
end


"""
    function get_geometry_center(cylinder::GadgetCylinder)

Returns the center of `cylinder`.
"""
get_geometry_center(cylinder::GadgetCylinder) = 1 // 2 .* (cylinder.pos_start .+ cylinder.pos_end)

"""
    get_geometry_box_corners(cylinder::GadgetCylinder)

Returns a tuple with the lower left and upper right corner of a box which contains the `cylinder`.
"""
function get_geometry_box_corners(cylinder::GadgetCylinder)
    min.(cylinder.pos_start, cylinder.pos_end) .- cylinder.radius,
    max.(cylinder.pos_start, cylinder.pos_end) .+ cylinder.radius
end


"""
    get_geometry_mask(cylinder::GadgetCylinder, pos::Matrix{T}) where T

Returns the indices of all particles contained in the `cylinder`.
"""
function get_geometry_mask(cylinder::GadgetCylinder, pos::Matrix{T}) where T

    # allocate mask array
    mask = Vector{Bool}(undef, size(pos,2))

    # check every position array if it's contained in the cylinder 
    # and update the mask array
    @inbounds for i = 1:size(pos,2)
        mask[i] = check_in_cylinder(pos[:,i], 
                                    cylinder.pos_start, 
                                    cylinder.pos_end, 
                                    cylinder.radius)
    end
   
    return findall(mask)
end
