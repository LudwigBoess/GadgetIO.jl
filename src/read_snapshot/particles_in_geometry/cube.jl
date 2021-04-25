"""
    struct GadgetCube{T} <: AbstractGadgetGeometry
        corner_lower_left::Vector{T}
        corner_upper_right::Vector{T}
    end

Defines a cube by `corner_lower_left` and `corner_upper_right`. 
To be used for [`read_particles_in_geometry`](@ref)
"""
struct GadgetCube{T} <: AbstractGadgetGeometry
    corner_lower_left::Vector{T}
    corner_upper_right::Vector{T}
end


"""
    get_geometry_box_corners(cube::GadgetCube)

Returns a `Tuple` with the lower left and upper right corner of a box which contains the `cube`.
"""
function get_geometry_box_corners(cube::GadgetCube)
    cylinder.pos_start .- cylinder.radius, cylinder.pos_end .+ cylinder.radius
end


"""
    get_geometry_mask(cube::GadgetCube, pos::Matrix{T}) where T

Returns the indices of all particles contained in the `cube`.
"""
function get_geometry_mask(cube::GadgetCube, pos::Matrix{T}) where T
    # by definition all particles are in the cube
    return collect(1:size(pos,2))
end