"""
    snap_to_dict(filename::String, try_info::Bool=true)

Reads whole snapshot into memory and returns a dictionary sorted by particle type
and block names. Be cautious with large snapshots!
"""
function snap_to_dict(filename::String, try_info::Bool=true)

        data = Dict()

        data["Header"] = head_to_dict(filename)

        h = head_to_obj(filename)

        if data["Header"]["snap_format"] == 2

            info = read_info(filename)
            blocks = print_blocks(filename, verbose=false)

            # delete "HEAD" and "INFO"
            delete_blocks = ["HEAD", "INFO"]

            for del_block ∈ delete_blocks
                deletepos = findfirst(blocks .== del_block)

                if !isnothing(deletepos)
                    deleteat!(blocks, deletepos)
                end
            end

            if typeof(info) == Array{InfoLine,1}

                block_positions = get_block_positions(filename)
                
                for parttype = 0:5
                    if data["Header"]["nall"][parttype+1] > 0
                        
                        # allocate dict for particle type
                        keyname = "PartType$parttype"
                        data[keyname] = Dict()

                        for block ∈ blocks
                            
                            block_info = info[getfield.(info, :block_name) .== block][1]

                            if block_info.is_present[parttype+1] == 1
                                data[keyname][block] = read_block(filename, block, 
                                                                  parttype=parttype,
                                                                  block_position=block_positions[block],
                                                                  info=block_info,
                                                                  h=h)
                            end

                        end
                    end
                end

            else

                error("No INFO block present! Snapshot can't be read in the lazy way.")

            end

        else
            error("Reading format 1 has been removed in version 0.3.0! If you need it please open an issue on GitHub!")
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
