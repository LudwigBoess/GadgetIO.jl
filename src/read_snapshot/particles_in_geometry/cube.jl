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
    cube.corner_lower_left, cube.corner_upper_right
end

"""
    function get_geometry_center(sphere::GadgetCube)

Returns the center of `cube`.
"""
get_geometry_center(cube::GadgetCube) =  1 // 2 .* cube.corner_lower_left .+ cube.corner_upper_right

"""
    get_geometry_mask(cube::GadgetCube, pos::Matrix{T}) where T

Returns the indices of all particles contained in the `cube`.
"""
function get_geometry_mask(cube::GadgetCube, pos::Matrix{T}) where T
    # by definition all particles are in the cube
    #return collect(1:size(pos,2))
    mask = @views @. ( ( cube.corner_lower_left[1] <= pos[1,:] <= cube.corner_upper_right[1] ) &
                       ( cube.corner_lower_left[2] <= pos[2,:] <= cube.corner_upper_right[2] ) &
                       ( cube.corner_lower_left[3] <= pos[3,:] <= cube.corner_upper_right[3] ) )

    return findall(mask)
end