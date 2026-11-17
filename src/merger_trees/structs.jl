"""
Merger-tree data structures.

Ports the on-disk layout used by the C codes L-BaseTree and L-HaloTrees
(Springel, Lemson). The binary layout of `MergerTreeHalo` matches the
`halo_data` C struct so that `read!(io, Vector{MergerTreeHalo})` decodes
existing `trees_NNN.X` files directly.
"""

"""
    const NO_HALO = Int32(-1)

Sentinel value used by pointer fields of [`MergerTreeHalo`](@ref)
(`Descendant`, `FirstProgenitor`, `NextProgenitor`,
`FirstHaloInFOFgroup`, `NextHaloInFOFgroup`) to indicate "no halo".
Matches the on-disk convention of the L-HaloTrees C code.
"""
const NO_HALO = Int32(-1)

"""
    struct MergerTreeHalo

A single node in a merger tree. Layout mirrors the `halo_data` C struct
used by `L-HaloTrees` (104 bytes, no padding):

| Field                  | Type             | Meaning                                                                |
| :--------------------- | :--------------- | :--------------------------------------------------------------------- |
| `Descendant`           | `Int32`          | in-tree index of descendant halo, `-1` if none                         |
| `FirstProgenitor`      | `Int32`          | in-tree index of largest progenitor, `-1` if none                      |
| `NextProgenitor`       | `Int32`          | next progenitor of the same descendant (linked list), `-1` if none     |
| `FirstHaloInFOFgroup`  | `Int32`          | index of central halo of this FOF group                                |
| `NextHaloInFOFgroup`   | `Int32`          | next satellite in the same FOF group, `-1` if none                     |
| `Len`                  | `Int32`          | number of particles in the (sub-)halo                                  |
| `M_Mean200`            | `Float32`        | mass enclosing mean overdensity 200 (only for centrals)                |
| `M_Crit200`            | `Float32`        | mass enclosing critical overdensity 200                                |
| `M_TopHat`             | `Float32`        | top-hat virial mass                                                    |
| `Pos`                  | `NTuple{3,F32}`  | position                                                               |
| `Vel`                  | `NTuple{3,F32}`  | velocity                                                               |
| `VelDisp`              | `Float32`        | velocity dispersion                                                    |
| `Vmax`                 | `Float32`        | maximum circular velocity                                              |
| `Spin`                 | `NTuple{3,F32}`  | spin vector                                                            |
| `MostBoundID`          | `Int64`          | particle ID of most bound particle                                     |
| `SnapNum`              | `Int32`          | snapshot in which the halo was found                                   |
| `FileNr`               | `Int32`          | subfind file containing this halo                                      |
| `SubhaloIndex`         | `Int32`          | subhalo index inside its subfind file                                  |
| `SubhalfMass`          | `Float32`        | half-mass radius / mass marker                                         |
"""
struct MergerTreeHalo
    Descendant::Int32
    FirstProgenitor::Int32
    NextProgenitor::Int32
    FirstHaloInFOFgroup::Int32
    NextHaloInFOFgroup::Int32
    Len::Int32
    M_Mean200::Float32
    M_Crit200::Float32
    M_TopHat::Float32
    Pos::NTuple{3,Float32}
    Vel::NTuple{3,Float32}
    VelDisp::Float32
    Vmax::Float32
    Spin::NTuple{3,Float32}
    MostBoundID::Int64
    SnapNum::Int32
    FileNr::Int32
    SubhaloIndex::Int32
    SubhalfMass::Float32
end

# sanity: keep this in sync with the C `halo_data` layout (104 bytes)
@assert isbitstype(MergerTreeHalo)
@assert sizeof(MergerTreeHalo) == 104

"""
    MergerTreeHalo()

Construct an empty (all-zero, dangling pointers set to `-1`) halo. Useful
when building trees from scratch.
"""
function MergerTreeHalo()
    return MergerTreeHalo(
        NO_HALO, NO_HALO, NO_HALO, NO_HALO, NO_HALO,
        Int32(0),
        0f0, 0f0, 0f0,
        (0f0, 0f0, 0f0), (0f0, 0f0, 0f0),
        0f0, 0f0,
        (0f0, 0f0, 0f0),
        Int64(0),
        Int32(0), Int32(0), Int32(0),
        0f0,
    )
end

"""
    struct MergerTree

One merger tree: a flat depth-first ordered vector of `MergerTreeHalo`.
All pointer fields (`Descendant`, `FirstProgenitor`, ...) are **0-based
indices into `halos`**, matching the on-disk convention. A value of `-1`
means "no halo".

Use `treeidx(tree, halo.Descendant)` to convert to a 1-based Julia index.
"""
struct MergerTree
    halos::Vector{MergerTreeHalo}
end

Base.length(t::MergerTree) = length(t.halos)
Base.getindex(t::MergerTree, i::Integer) = t.halos[i]
Base.iterate(t::MergerTree, state...) = iterate(t.halos, state...)
Base.firstindex(t::MergerTree) = firstindex(t.halos)
Base.lastindex(t::MergerTree) = lastindex(t.halos)

"""
    treeidx(tree, ptr)

Convert a 0-based pointer field (`-1` = none) to a 1-based Julia index
into `tree.halos`, or `nothing` if the pointer is `-1`.
"""
@inline function treeidx(::MergerTree, ptr::Integer)
    return ptr < 0 ? nothing : Int(ptr) + 1
end

"""
    struct MergerTreeFile

The contents of a single `trees_NNN.X` file: a collection of independent
merger trees.

| Field          | Type                 | Meaning                                  |
| :------------- | :------------------- | :--------------------------------------- |
| `trees`        | `Vector{MergerTree}` | one entry per tree                       |
| `ntrees`       | `Int32`              | number of trees                          |
| `tot_nhalos`   | `Int32`              | total halos summed over all trees        |
"""
struct MergerTreeFile
    trees::Vector{MergerTree}
    ntrees::Int32
    tot_nhalos::Int32
end

Base.length(f::MergerTreeFile) = length(f.trees)
Base.getindex(f::MergerTreeFile, i::Integer) = f.trees[i]
Base.iterate(f::MergerTreeFile, state...) = iterate(f.trees, state...)
