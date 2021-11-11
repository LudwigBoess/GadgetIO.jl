using LinearAlgebra

"""
    filter_cube(snap_file::String, corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real}, 
                          parttype::Integer)

Reads positions from `snap_file` and returns the indices of particles contained in a box defined by `corner_lowerleft` and `corner_upperright`.
"""
function filter_cube(snap_file::String, corner_lowerleft::Array{<:Real}, corner_upperright::Array{<:Real}; parttype::Integer=0)

    # parttype not in subfile
    if check_info(snap_file, "POS").is_present[parttype+1] == 0
        return Int[]
    end

    # read positions from file
    pos = read_snap(snap_file, "POS", parttype)

    # double check if the corners are set correctly
    cube_lower = Array{eltype(corner_lowerleft[1]),1}(undef, 3)
    cube_upper = Array{eltype(corner_lowerleft[1]),1}(undef, 3)
    @inbounds for dim = 1:3
        cube_lower[dim] = ((corner_lowerleft[dim] < corner_upperright[dim]) ? corner_lowerleft[dim] : corner_upperright[dim])
        cube_upper[dim] = ((corner_lowerleft[dim] > corner_upperright[dim]) ? corner_lowerleft[dim] : corner_upperright[dim])
    end

    return findall( @views ( ( cube_lower[1] .<= pos[1,:] .<= cube_upper[1] ) .& 
                             ( cube_lower[2] .<= pos[2,:] .<= cube_upper[2] ) .&
                             ( cube_lower[3] .<= pos[3,:] .<= cube_upper[3] )) )
end

"""
    check_in_cylinder(x, pt1, pt2, r)

Checks if a 3D point `x` is in a cylinder defined by its endpoints `pt1` and `pt2` and radius `r`.
"""
function check_in_cylinder(x, pt1, pt2, r)
    # https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder
    ( ( x .- pt1 ) ⋅ ( pt2 .- pt1 ) ≥ 0.0                             ) & 
    ( ( x .- pt2 ) ⋅ ( pt2 .- pt1 ) ≤ 0.0                             ) & 
    ( ( norm( (x .- pt1) × (pt2 .- pt1 ) ) / norm( pt2 .- pt1 ) ) ≤ r )
end

"""
    filter_cylinder(filename::String, pt1::Array{<:Real}, pt2::Array{<:Real}, r::Real)

Reads the positions contained in a file and returns the indices of particles contained in a cylinder defined by the endpoints `pt1` and `pt2` and radius `r`.
"""
function filter_cylinder(filename::String, pt1::Array{<:Real}, pt2::Array{<:Real}, r::Real; parttype::Integer=0)

    in_cube = filter_cube(filename, pt1, pt2, parttype=parttype )

    if size(in_cube,1) > 0
        check_in_cylinder_helper(x) = check_in_cylinder(x, pt1, pt2, r)
        # read positions from file
        pos = read_snap(filename, "POS", parttype)

        pos_cylinder = @views pos[:, in_cube]
        is_in_cylinder = mapslices(check_in_cylinder_helper, pos_cylinder, dims=1)
        return in_cube[findall( is_in_cylinder .== true )]
    else
        return in_cube
    end
end



"""
    filter_sphere(filename::String, center::Array{<:Real}, r::Real; parttype::Integer=0)

Reads the positions contained in a file and returns the indices of particles contained in a cylinder defined by the endpoints `pt1` and `pt2` and radius `r`.
"""
function filter_sphere(filename::String, center::Array{<:Real}, r::Real; parttype::Integer=0)

    in_cube = filter_cube(filename, center .- r, center .+ r, parttype=parttype)

    if size(in_cube,1) > 0

        distance = @views @. √( (pos[1, in_cube] - center[1])^2 +
                                (pos[2, in_cube] - center[2])^2 + 
                                (pos[3, in_cube] - center[3])^2 )

        return in_cube[findall( distance .<= r )]
    else
        return in_cube
    end
end
