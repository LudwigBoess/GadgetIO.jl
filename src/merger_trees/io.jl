"""
Binary I/O for L-BaseTree and L-HaloTrees output files.

`trees_NNN.X` layout (L-HaloTrees):
```
Int32  ntrees
Int32  tot_nhalos
Int32  nhalos_per_tree[ntrees]
MergerTreeHalo halos[tot_nhalos]   # depth-first, concatenated
```

`sub_desc_NNN` layout (L-BaseTree):
```
Int32  ntot                                  # number of subhalos in source snap
Int32  descendant_haloindex[ntot]            # -1 if no descendant
Int32  descendant_snapnum[ntot]              # -1 if no descendant
[Int64 most_bound_id[ntot]]                  # only with LGADGET3
```
"""

"""
    read_merger_tree_file(filename) -> MergerTreeFile

Read a single `trees_NNN.X` file written by L-HaloTrees. All pointer
fields stay 0-based to match the original convention.
"""
function read_merger_tree_file(filename::AbstractString)
    open(filename, "r") do io
        ntrees = read(io, Int32)
        tot_nhalos = read(io, Int32)

        nhalos_per_tree = Vector{Int32}(undef, ntrees)
        read!(io, nhalos_per_tree)

        trees = Vector{MergerTree}(undef, ntrees)
        for t in 1:ntrees
            n = Int(nhalos_per_tree[t])
            halos = Vector{MergerTreeHalo}(undef, n)
            read!(io, halos)
            trees[t] = MergerTree(halos)
        end

        return MergerTreeFile(trees, ntrees, tot_nhalos)
    end
end

"""
    read_merger_trees(treedir; snapnum, files=:all) -> Vector{MergerTreeFile}

Read every `trees_<snapnum>.X` file in `treedir`. Pass `files=0:N-1` to
load a subset.
"""
function read_merger_trees(treedir::AbstractString;
                           snapnum::Integer,
                           files = :all)
    base = joinpath(treedir, @sprintf("trees_%03d", snapnum))

    if files === :all
        # discover files of the form trees_NNN.X
        filenames = sort!(filter(f -> startswith(f, basename(base) * "."),
                                 readdir(treedir)))
        files = [parse(Int, split(f, '.')[end]) for f in filenames]
    end

    return [read_merger_tree_file(string(base, '.', i)) for i in files]
end

"""
    write_merger_tree_file(filename, file::MergerTreeFile)

Write a `MergerTreeFile` back to disk in the L-HaloTrees binary format.
"""
function write_merger_tree_file(filename::AbstractString, file::MergerTreeFile)
    open(filename, "w") do io
        write(io, Int32(file.ntrees))
        write(io, Int32(file.tot_nhalos))

        nhalos_per_tree = Int32[length(t) for t in file.trees]
        write(io, nhalos_per_tree)

        for t in file.trees
            write(io, t.halos)
        end
    end
    return filename
end


"""
    struct SubDesc

Per-source-snapshot descendant information written by L-BaseTree.

| Field                    | Meaning                                                |
| :----------------------- | :----------------------------------------------------- |
| `nsubhalos`              | number of subhalos in the source snapshot              |
| `descendant_haloindex`   | descendant index in the target snap (0-based, -1 none) |
| `descendant_snapnum`     | snapshot number of the descendant, -1 if none          |
| `most_bound_id`          | optional, LGADGET3 only                                |
"""
struct SubDesc
    nsubhalos::Int32
    descendant_haloindex::Vector{Int32}
    descendant_snapnum::Vector{Int32}
    most_bound_id::Vector{Int64}
end

"""
    read_sub_desc(filename; lgadget3=false) -> SubDesc

Read a `sub_desc_NNN` file written by L-BaseTree. Set `lgadget3=true`
when the file was produced with the LGADGET3 most-bound-ID block.
"""
function read_sub_desc(filename::AbstractString; lgadget3::Bool=false)
    open(filename, "r") do io
        n = read(io, Int32)
        hidx = Vector{Int32}(undef, n)
        read!(io, hidx)
        hsnap = Vector{Int32}(undef, n)
        read!(io, hsnap)

        mbid = Int64[]
        if lgadget3
            mbid = Vector{Int64}(undef, n)
            read!(io, mbid)
        end

        return SubDesc(n, hidx, hsnap, mbid)
    end
end

"""
    write_sub_desc(filename, d::SubDesc; lgadget3=false)

Write a `sub_desc_NNN` file in the L-BaseTree binary format.
"""
function write_sub_desc(filename::AbstractString, d::SubDesc; lgadget3::Bool=false)
    open(filename, "w") do io
        write(io, Int32(d.nsubhalos))
        write(io, d.descendant_haloindex)
        write(io, d.descendant_snapnum)
        if lgadget3
            @assert length(d.most_bound_id) == d.nsubhalos
            write(io, d.most_bound_id)
        end
    end
    return filename
end
