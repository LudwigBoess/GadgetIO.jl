"""
    snap_to_dict(filename::String, try_info::Bool=true)

Reads whole snapshot into memory and returns a dictionary sorted by particle type
and block names. Be cautious with large snapshots!
"""
function snap_to_dict(filename::String, try_info::Bool=true)

        data = Dict()

        data["Header"] = head_to_dict(filename)

        if data["Header"]["snap_format"] == 2

            if try_info

                    info = read_info(filename)

                    if typeof(info) == Array{InfoLine,1}

                        data = snap_2_d_info(filename, data, info)

                    else

                        data = snap_2_d(filename, data)

                    end

            else

                data = snap_2_d(filename, data)

            end
        else

            data = snap_1_d(filename, data)

        end

        return data
end



"""
    read_snap(filename::String [, blockname::String="", parttype::Integer=-1] )

Wrapper function to read snapshot in various ways:
filename only: returns the entire snapshot as a dictionary.
blockname: Returns only that block. If parttype specified only for that
particle type.

# Examples
```jldoctest
julia> gas_pos = read_block(filename, "POS", 0)
[...]
```

"""
function read_snap(filename::String, blockname::String="", parttype::Integer=-1)

    # default: return entire snapshot as dictionary
    if (parttype == -1) & (blockname == "")
        return snap_to_dict(filename)
    end

    if blockname != ""
        return read_block(filename, blockname, parttype=parttype)
    end

end
