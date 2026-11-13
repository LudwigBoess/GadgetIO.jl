"""
Reader for **consistent-trees** (Behroozi et al.) merger-tree output.

Consistent-trees writes ASCII `tree_*.dat` files inside the
`<rockstar_run>/trees/` directory. Each file holds many independent
merger trees in the layout

```
#scale(0) id(1) desc_scale(2) desc_id(3) ... Tidal_ID(36) <extra columns>
#Omega_M = ...; Omega_L = ...; h0 = ...
#  ... more comment lines describing every column ...
<ntrees>                         # single integer line, number of trees in file
#tree <root_id>                  # delimiter, one per tree
<data row>                       # depth-first ordered halos of that tree
<data row>
...
#tree <root_id>
...
```

The first non-comment line is a bare integer giving the number of trees;
each subsequent tree is introduced by a `#tree <root_id>` comment line and
followed by its halos in **depth-first order** (the root, whose
`desc_id == -1`, comes first).

This reader converts each consistent-trees tree into a
[`MergerTree`](@ref) that is bit-for-bit compatible with the
L-HaloTrees-based trees produced elsewhere in this package, so that the
same downstream tooling (`root_index`, `walk_main_branch`,
`num_progenitors`, and `MakieHelper.plot_merger_tree`) works unchanged.

Column positions are **not** hardcoded: the leading `#scale(0) id(1) ...`
header line is parsed into a `name → 0-based column index` map, so the
reader tolerates a different `MASS_DEFINITION` (e.g. `"m200b"` instead of
`"vir"`) or a different set of trailing extra parameters. Optional columns
(`M200b`, `M200c`, ...) are used when present and gracefully fall back to
`Mvir` otherwise.

# CT column → `MergerTreeHalo` mapping

| `MergerTreeHalo` field | consistent-trees source                                    |
| :--------------------- | :--------------------------------------------------------- |
| `Descendant`           | `desc_id` mapped to local 0-based index (`-1` = root)      |
| `FirstProgenitor`      | built by [`set_progenitor_pointers!`](@ref) from `Descendant` |
| `NextProgenitor`       | built by [`set_progenitor_pointers!`](@ref)                |
| `FirstHaloInFOFgroup`  | central via `upid` (`upid==-1` ⇒ self), else `NO_HALO`     |
| `NextHaloInFOFgroup`   | `NO_HALO` (not reconstructed; unused by plotting)          |
| `Len`                  | `round(Int32, Mvir / particle_mass)`                       |
| `M_Mean200`            | `M200b` if present, else `Mvir`                            |
| `M_Crit200`            | `M200c` if present, else `Mvir`                            |
| `M_TopHat`             | `Mvir`                                                      |
| `Pos`                  | `(x, y, z)`                                                 |
| `Vel`                  | `(vx, vy, vz)`                                              |
| `VelDisp`              | `vrms`                                                      |
| `Vmax`                 | `vmax`                                                      |
| `Spin`                 | `(Jx, Jy, Jz)`                                              |
| `MostBoundID`          | `id`                                                        |
| `SnapNum`              | `Snap_idx`                                                  |
| `FileNr`               | `0`                                                         |
| `SubhaloIndex`         | `Orig_halo_ID`                                              |
| `SubhalfMass`          | `rs`                                                        |

Masses are in `Msun/h`, positions in `Mpc/h` (comoving), exactly as
consistent-trees emits them.
"""

# default high-resolution dark-matter particle mass (Msun/h) used to turn
# Mvir into an integer particle count for MakieHelper's `log10(Len)` default.
const CT_DEFAULT_PARTICLE_MASS = 3.873031e7

"""
    parse_ct_header(line) -> Dict{String,Int}

Parse a consistent-trees column header line of the form

```
#scale(0) id(1) desc_scale(2) desc_id(3) num_prog(4) ... Tidal_ID(36) Rs_Klypin(37) ...
```

into a map of (lower-cased) column **name → 0-based column index**. Each
token looks like `name(idx)`; the parenthesised index is authoritative
(so the map is robust even if tokens are reordered). Tokens without an
explicit `(idx)` are assigned their running position as a fallback.

Punctuation that only decorates a name in the header but never appears in
the lookup keys we query (`?`, `/`, `|`, `[`, `]`, ...) is preserved in
the returned keys verbatim; query through [`ct_col`](@ref), which
lower-cases and tries a couple of spellings.
"""
function parse_ct_header(line::AbstractString)
    cols = Dict{String,Int}()
    s = lstrip(line)
    startswith(s, "#") && (s = lstrip(s[2:end]))
    running = 0
    for tok in split(s)
        isempty(tok) && continue
        name = tok
        idx  = running
        m = match(r"^(.*)\((\d+)\)$", tok)
        if m !== nothing
            name = m.captures[1]
            idx  = parse(Int, m.captures[2])
        end
        cols[lowercase(name)] = idx
        running += 1
    end
    return cols
end

"""
    ct_col(cols, names...) -> Union{Int,Nothing}

Look up the 0-based column index for the first of `names` that is present
in the `name → index` map produced by [`parse_ct_header`](@ref).
Comparison is case-insensitive. Returns `nothing` if none match, which
callers use to fall back on optional columns.
"""
function ct_col(cols::AbstractDict{String,Int}, names::AbstractString...)
    for nm in names
        idx = get(cols, lowercase(nm), nothing)
        idx === nothing || return idx
    end
    return nothing
end

# fetch a numeric field from a split data row by 0-based column index,
# with a default for absent/optional columns.
@inline function _ct_get(fields::AbstractVector{<:AbstractString},
                         idx::Union{Int,Nothing}, ::Type{T}, default::T) where {T}
    idx === nothing && return default
    (idx + 1) > length(fields) && return default
    return parse(T, fields[idx + 1])
end

"""
    _build_ct_tree(rows, cols; particle_mass) -> MergerTree

Convert the depth-first ordered consistent-trees data `rows` of a single
tree (each `rows[k]` already split into whitespace-separated fields) into
a [`MergerTree`](@ref). Builds the `ct_id → local 0-based index` map,
fills every `MergerTreeHalo` field from the column map `cols`, then calls
[`set_progenitor_pointers!`](@ref) to construct the
`FirstProgenitor`/`NextProgenitor` linked lists exactly as the
L-HaloTrees path does.
"""
function _build_ct_tree(rows::AbstractVector{<:AbstractVector{<:AbstractString}},
                        cols::AbstractDict{String,Int};
                        particle_mass::Real = CT_DEFAULT_PARTICLE_MASS)

    n = length(rows)

    # resolve column indices once for the whole tree
    c_id       = ct_col(cols, "id")
    c_descid   = ct_col(cols, "desc_id", "descid")
    c_upid     = ct_col(cols, "upid")
    c_mvir     = ct_col(cols, "mvir")
    c_m200b    = ct_col(cols, "m200b")
    c_m200c    = ct_col(cols, "m200c")
    c_rs       = ct_col(cols, "rs")
    c_vrms     = ct_col(cols, "vrms")
    c_vmax     = ct_col(cols, "vmax")
    c_x        = ct_col(cols, "x")
    c_y        = ct_col(cols, "y")
    c_z        = ct_col(cols, "z")
    c_vx       = ct_col(cols, "vx")
    c_vy       = ct_col(cols, "vy")
    c_vz       = ct_col(cols, "vz")
    c_jx       = ct_col(cols, "jx")
    c_jy       = ct_col(cols, "jy")
    c_jz       = ct_col(cols, "jz")
    c_snap     = ct_col(cols, "snap_idx", "snap_num")
    c_orig     = ct_col(cols, "orig_halo_id")

    c_id   === nothing && error("consistent-trees header is missing the `id` column")
    c_descid === nothing && error("consistent-trees header is missing the `desc_id` column")
    c_mvir === nothing && error("consistent-trees header is missing the `Mvir` column")

    # 1) ct_id -> local 0-based index, preserving the CT depth-first order.
    id_to_local = Dict{Int64,Int}()
    sizehint!(id_to_local, n)
    @inbounds for k in 1:n
        ctid = parse(Int64, rows[k][c_id + 1])
        id_to_local[ctid] = k - 1
    end

    halos = Vector{MergerTreeHalo}(undef, n)

    @inbounds for k in 1:n
        f = rows[k]

        ctid   = parse(Int64, f[c_id + 1])
        descid = _ct_get(f, c_descid, Int64, Int64(-1))

        # Descendant: map CT desc_id to a *local* 0-based index. Root has
        # desc_id == -1 -> NO_HALO. A desc_id that points outside this
        # tree (should not happen for well-formed CT output) is treated as
        # a root for safety.
        descendant = NO_HALO
        if descid != -1
            li = get(id_to_local, descid, nothing)
            descendant = li === nothing ? NO_HALO : Int32(li)
        end

        mvir  = _ct_get(f, c_mvir,  Float32, 0f0)
        m200b = c_m200b === nothing ? mvir : _ct_get(f, c_m200b, Float32, mvir)
        m200c = c_m200c === nothing ? mvir : _ct_get(f, c_m200c, Float32, mvir)

        len = round(Int32, Float64(mvir) / Float64(particle_mass))

        x  = _ct_get(f, c_x,  Float32, 0f0)
        y  = _ct_get(f, c_y,  Float32, 0f0)
        z  = _ct_get(f, c_z,  Float32, 0f0)
        vx = _ct_get(f, c_vx, Float32, 0f0)
        vy = _ct_get(f, c_vy, Float32, 0f0)
        vz = _ct_get(f, c_vz, Float32, 0f0)
        jx = _ct_get(f, c_jx, Float32, 0f0)
        jy = _ct_get(f, c_jy, Float32, 0f0)
        jz = _ct_get(f, c_jz, Float32, 0f0)

        vrms = _ct_get(f, c_vrms, Float32, 0f0)
        vmax = _ct_get(f, c_vmax, Float32, 0f0)
        rs   = _ct_get(f, c_rs,   Float32, 0f0)

        snap = _ct_get(f, c_snap, Int32, Int32(0))

        # Orig_halo_ID is an int64 in CT and can exceed Int32 range in big
        # boxes; parse wide and wrap into the Int32 SubhaloIndex slot rather
        # than throwing on overflow (the field is informational for plotting).
        orig64 = _ct_get(f, c_orig, Int64, Int64(0))
        orig   = orig64 % Int32

        # FOF linkage: best-effort central via `upid`. upid==-1 marks a
        # distinct (central) halo -> point FirstHaloInFOFgroup at itself.
        # Satellites point at their central's local index when it is in
        # this tree. NextHaloInFOFgroup is left as NO_HALO; neither field
        # is consulted by MakieHelper's layout, so this is purely for
        # bookkeeping fidelity.
        first_in_fof = NO_HALO
        upid = _ct_get(f, c_upid, Int64, Int64(-1))
        if upid == -1
            first_in_fof = Int32(k - 1)
        else
            li = get(id_to_local, upid, nothing)
            first_in_fof = li === nothing ? NO_HALO : Int32(li)
        end

        halos[k] = MergerTreeHalo(
            descendant,
            NO_HALO,          # FirstProgenitor — filled by set_progenitor_pointers!
            NO_HALO,          # NextProgenitor  — filled by set_progenitor_pointers!
            first_in_fof,
            NO_HALO,          # NextHaloInFOFgroup
            len,
            m200b,            # M_Mean200
            m200c,            # M_Crit200
            mvir,             # M_TopHat
            (x, y, z),
            (vx, vy, vz),
            vrms,
            vmax,
            (jx, jy, jz),
            ctid,             # MostBoundID — CT unique halo id
            snap,
            Int32(0),         # FileNr
            orig,             # SubhaloIndex — Orig_halo_ID
            rs,               # SubhalfMass — scale radius marker
        )
    end

    # Build FirstProgenitor / NextProgenitor exactly like the L-HaloTrees
    # path: MMP (largest Len) sits at the head of each descendant's chain.
    # Because Len ∝ Mvir here, this reproduces CT's mass ordering and the
    # `mmp?` flag, which is what is_mmp / num_progenitors / tree_layout
    # expect.
    set_progenitor_pointers!(halos)

    return MergerTree(halos)
end


"""
    read_consistent_trees(filename; particle_mass=CT_DEFAULT_PARTICLE_MASS)
        -> Vector{MergerTree}

Read a single consistent-trees `tree_*.dat` file and return one
[`MergerTree`](@ref) per tree it contains, ready for analysis or for
`MakieHelper.plot_merger_tree`.

The column layout is taken from the `#scale(0) id(1) ...` header line
(see [`parse_ct_header`](@ref)), so a non-`vir` `MASS_DEFINITION` or extra
trailing parameters are handled automatically. `particle_mass` (in
`Msun/h`) sets the conversion `Len = round(Mvir / particle_mass)` so that
MakieHelper's default `mass_of = h -> log10(Len)` is sensible; for science
plots prefer passing an explicit mass accessor, e.g.
`mass_of = h -> log10(h.M_Crit200)`.

# Example
```julia
using GadgetIO, MakieHelper, CairoMakie

trees = read_consistent_trees("rockstar_run/trees/tree_0_0_0.dat")
@show length(trees)

# real scale factors for the x-axis (see scale_factor_map)
scaledict = scale_factor_map("rockstar_run/trees/tree_0_0_0.dat")

fig, ax = plot_merger_tree(trees[1];
                           mass_of = h -> log10(h.M_Crit200),
                           snap_to_scale = s -> scaledict[s])
```
"""
function read_consistent_trees(filename::AbstractString;
                               particle_mass::Real = CT_DEFAULT_PARTICLE_MASS)

    cols::Union{Dict{String,Int},Nothing} = nothing
    trees = MergerTree[]

    # per-tree accumulator of split data rows; reused (emptied) per tree so
    # we never reassign a captured local across the file/closure boundary.
    current = Vector{Vector{SubString{String}}}()
    in_tree = false

    # flush the rows accumulated for the current tree into `trees`
    flush_tree!() =
        if in_tree && !isempty(current)
            push!(trees, _build_ct_tree(current, cols; particle_mass))
            empty!(current)
        end

    open(filename, "r") do io
        for line in eachline(io)
            s = lstrip(line)
            isempty(s) && continue

            if startswith(s, "#")
                # the first header token starting with "#scale" carries the
                # column layout; "#tree <id>" (lower-case, space, then the root
                # id) starts a new tree.  Match it case-sensitively and require
                # the trailing space so the description line "#Tree_root_ID: ..."
                # (capital T, underscore) is NOT mistaken for a tree delimiter —
                # that false match would otherwise swallow the bare ntrees count
                # line into a bogus one-column tree.
                if cols === nothing && startswith(lowercase(s), "#scale")
                    cols = parse_ct_header(s)
                elseif startswith(s, "#tree ")
                    flush_tree!()
                    in_tree = true
                end
                continue
            end

            # non-comment line: either the bare ntrees count (before the
            # first "#tree") or a halo data row.
            if !in_tree
                # ntrees count line — informational only; we count trees
                # from the "#tree" delimiters instead, so just skip it.
                continue
            end

            cols === nothing &&
                error("consistent-trees file $(filename) has data before its `#scale(...)` header")

            push!(current, split(s))
        end
    end
    flush_tree!()

    return trees
end


"""
    read_consistent_trees_dir(treedir; particle_mass=CT_DEFAULT_PARTICLE_MASS)
        -> Vector{MergerTree}

Read every `tree_*.dat` file in a consistent-trees `trees/` directory and
return the concatenation of all trees they contain (one
[`MergerTree`](@ref) per tree). Files are processed in sorted filename
order. See [`read_consistent_trees`](@ref) for the per-file reader and the
column mapping.

# Example
```julia
using GadgetIO
trees = read_consistent_trees_dir("rockstar_run/trees")
@show length(trees)
```
"""
function read_consistent_trees_dir(treedir::AbstractString;
                                   particle_mass::Real = CT_DEFAULT_PARTICLE_MASS)
    files = sort!(filter(f -> startswith(f, "tree_") && endswith(f, ".dat"),
                         readdir(treedir)))
    isempty(files) &&
        error("no tree_*.dat files found in $(treedir)")

    trees = MergerTree[]
    for f in files
        append!(trees, read_consistent_trees(joinpath(treedir, f);
                                             particle_mass))
    end
    return trees
end


"""
    read_most_massive_tree(filename; mass_field="Mvir",
                           particle_mass=CT_DEFAULT_PARTICLE_MASS) -> MergerTree

Read only the single merger tree whose **root** halo (the present-day halo, the
one with `desc_id == -1`) has the largest `mass_field`, without building a
`MergerTree` for every other tree in the file.

[`read_consistent_trees`](@ref) materialises one `MergerTree` per tree in the
file (>10^5 for a full box, ~1 GB). When you only want the assembly history of
the most massive halo, this streams the file once, keeps just the rows of the
current best candidate, and converts that single tree at the end — a small
fraction of the time and memory.

"Most massive" is judged by the **root** halo's `mass_field` (default `"Mvir"`,
i.e. the present-day virial mass), matched case-insensitively against the file's
column header (so `"m200c"`, `"m200b"`, ... also work). consistent-trees emits
each tree's root first, so only the root rows are parsed during the scan.

# Example
```julia
using GadgetIO, MakieHelper, CairoMakie
tf   = "rockstar_run/trees/tree_0_0_0.dat"
tree = read_most_massive_tree(tf)
amap = scale_factor_map(tf)
fig, ax = plot_merger_tree(tree; mass_of = h -> log10(h.M_Crit200),
                                 snap_to_scale = s -> amap[s])
```

To instead pick the most massive halo across a whole `trees/` directory, take
the `argmax` of `read_most_massive_tree` over the `tree_*.dat` files.
"""
function read_most_massive_tree(filename::AbstractString;
                                mass_field::AbstractString = "Mvir",
                                particle_mass::Real = CT_DEFAULT_PARTICLE_MASS)
    cols::Union{Dict{String,Int},Nothing} = nothing
    c_mass::Union{Int,Nothing}   = nothing
    c_descid::Union{Int,Nothing} = nothing

    best_lines = String[]      # raw rows of the best tree found so far
    best_mass  = -Inf
    current    = String[]      # raw rows of the tree currently being read
    in_tree    = false

    # Root mass of the tree buffered in `lines`: the row with desc_id == -1.
    # consistent-trees emits the root first, so check the first row, with a
    # scan as a safety net if a variant ever reorders rows.
    root_mass(lines) = begin
        c_mass === nothing && error("column `$mass_field` not found in $(filename)")
        f = split(first(lines))
        (c_descid === nothing || (length(f) > c_descid && f[c_descid + 1] == "-1")) &&
            return parse(Float64, f[c_mass + 1])
        for l in lines
            g = split(l)
            length(g) > c_descid && g[c_descid + 1] == "-1" &&
                return parse(Float64, g[c_mass + 1])
        end
        return parse(Float64, f[c_mass + 1])
    end

    finish_tree!() =
        if !isempty(current)
            m = root_mass(current)
            if m > best_mass
                best_mass  = m
                best_lines = current     # keep this tree's rows ...
                current    = String[]    # ... and start a fresh buffer
            else
                empty!(current)          # discard, reuse buffer
            end
        end

    open(filename, "r") do io
        for line in eachline(io)
            s = lstrip(line)
            isempty(s) && continue
            if startswith(s, "#")
                if cols === nothing && startswith(lowercase(s), "#scale")
                    cols     = parse_ct_header(s)
                    c_mass   = ct_col(cols, mass_field)
                    c_descid = ct_col(cols, "desc_id", "descid")
                elseif startswith(s, "#tree ")
                    finish_tree!()
                    in_tree = true
                end
                continue
            end
            in_tree || continue
            cols === nothing &&
                error("consistent-trees file $(filename) has data before its `#scale(...)` header")
            push!(current, s)
        end
    end
    finish_tree!()

    isempty(best_lines) && error("no trees found in $(filename)")
    rows = [split(l) for l in best_lines]
    return _build_ct_tree(rows, cols; particle_mass)
end


"""
    scale_factor_map(filename) -> Dict{Int,Float64}

Build a `Snap_idx → scale` map from a consistent-trees `tree_*.dat` file
by reading the `scale` (column 0) and `Snap_idx` (column 31) values of
every data row. consistent-trees stores a real scale factor per halo, so
this lets you put the merger-tree plot on a physical scale-factor axis
instead of MakieHelper's `(snap+1)/(last_snap+1)` heuristic.

Pass the result to `plot_merger_tree` via the `snap_to_scale` keyword:

```julia
scaledict = scale_factor_map("rockstar_run/trees/tree_0_0_0.dat")
plot_merger_tree(trees[1];
                 snap_to_scale = s -> scaledict[s],
                 mass_of       = h -> log10(h.M_Crit200))
```

If a snapshot index appears with (slightly) differing scale factors across
halos, the last one read wins; in practice all halos of a given snapshot
share the same scale.
"""
function scale_factor_map(filename::AbstractString)
    scaledict = Dict{Int,Float64}()
    cols::Union{Dict{String,Int},Nothing} = nothing

    open(filename, "r") do io
        for line in eachline(io)
            s = lstrip(line)
            isempty(s) && continue

            if startswith(s, "#")
                if cols === nothing && startswith(lowercase(s), "#scale")
                    cols = parse_ct_header(s)
                end
                continue
            end

            cols === nothing && continue  # skip the ntrees count line
            fields = split(s)
            c_scale = ct_col(cols, "scale")
            c_snap  = ct_col(cols, "snap_idx", "snap_num")
            (c_scale === nothing || c_snap === nothing) && continue
            (length(fields) <= max(c_scale, c_snap)) && continue

            snap  = parse(Int, fields[c_snap + 1])
            scale = parse(Float64, fields[c_scale + 1])
            scaledict[snap] = scale
        end
    end

    return scaledict
end


# ===========================================================================
# Forests
# ===========================================================================
#
# consistent-trees groups its trees into *forests*: sets of trees whose halos
# are dynamically associated (they exchange particles / pass through each
# other) but whose z=0 roots are still distinct halos. An ongoing,
# not-yet-coalesced major merger therefore appears as *two roots in the same
# forest* — the main cluster and the infalling companion — rather than as a
# single merger node inside one tree. Plotting a single `MergerTree` can never
# show such a merger; the whole forest must be read.
#
# `forests.list`  : `#TreeRootID ForestID`           (one row per tree root)
# `locations.dat` : `#TreeRootID FileID Offset Filename Num_Halos Num_Bytes`
#                   `Offset` is the byte offset of the tree's first *data row*
#                   (immediately after its `#tree <id>` delimiter line).

"""
    CTTreeLocation

Where a single consistent-trees tree lives on disk, as recorded in
`locations.dat`: which `Filename` (relative to the `trees/` dir), the byte
`offset` of its first data row, the halo count `nhalos` and byte length
`nbytes`.
"""
struct CTTreeLocation
    fileid::Int
    offset::Int
    filename::String
    nhalos::Int
    nbytes::Int
end

"""
    read_forests_list(path) -> Dict{Int64,Int64}

Parse a consistent-trees `forests.list` into a `TreeRootID => ForestID` map.
"""
function read_forests_list(path::AbstractString)
    d = Dict{Int64,Int64}()
    open(path, "r") do io
        for line in eachline(io)
            s = lstrip(line)
            (isempty(s) || startswith(s, "#")) && continue
            f = split(s)
            length(f) < 2 && continue
            d[parse(Int64, f[1])] = parse(Int64, f[2])
        end
    end
    return d
end

"""
    read_tree_locations(path) -> Dict{Int64,CTTreeLocation}

Parse a consistent-trees `locations.dat` into a `TreeRootID => CTTreeLocation`
map, used to seek directly to an individual tree inside a (possibly large)
`tree_*.dat` file instead of streaming the whole file.
"""
function read_tree_locations(path::AbstractString)
    d = Dict{Int64,CTTreeLocation}()
    open(path, "r") do io
        for line in eachline(io)
            s = lstrip(line)
            (isempty(s) || startswith(s, "#")) && continue
            f = split(s)
            length(f) < 6 && continue
            d[parse(Int64, f[1])] = CTTreeLocation(parse(Int, f[2]),
                                                   parse(Int, f[3]),
                                                   String(f[4]),
                                                   parse(Int, f[5]),
                                                   parse(Int, f[6]))
        end
    end
    return d
end

# read the `#scale(0) ...` column header of a CT tree file once.
# NB: do not `return` from inside the `open(...) do ... end` block — that
# returns from the closure, not from `_read_ct_header`; capture and return after.
function _read_ct_header(filename::AbstractString)
    cols = nothing
    open(filename, "r") do io
        for line in eachline(io)
            s = lstrip(line)
            if startswith(lowercase(s), "#scale")
                cols = parse_ct_header(s)
                break
            end
        end
    end
    cols === nothing && error("no `#scale(...)` header line found in $filename")
    return cols
end

# read one tree whose first data row starts at byte `offset` (the
# `locations.dat` convention). Robust to the offset landing either on the
# first data row or on the `#tree <id>` delimiter line: skip a leading
# delimiter, collect data rows, and stop at the *next* `#tree`.
function _read_ct_tree_at(io::IO, offset::Integer,
                          cols::AbstractDict{String,Int};
                          particle_mass::Real)
    seek(io, offset)
    rows = Vector{Vector{SubString{String}}}()
    while !eof(io)
        s = lstrip(readline(io))
        isempty(s) && continue
        if startswith(s, "#tree ")
            isempty(rows) && continue   # leading delimiter of *this* tree
            break                        # delimiter of the *next* tree
        end
        startswith(s, "#") && continue
        push!(rows, split(s))
    end
    return _build_ct_tree(rows, cols; particle_mass)
end

"""
    struct MergerForest

A consistent-trees *forest*: the set of `MergerTree`s whose roots share a
`ForestID`. Iterating / indexing a `MergerForest` yields its `MergerTree`s
(largest first). `root_ids[k]` is the consistent-trees `TreeRootID` of
`trees[k]`.

Use [`build_merger_forest`](@ref) to construct one and
`MakieHelper.plot_merger_forest` to render every tree in one figure, which is
what reveals an ongoing (not-yet-coalesced) merger as two converging branches.
"""
struct MergerForest
    forest_id::Int64
    trees::Vector{MergerTree}
    root_ids::Vector{Int64}
end

Base.length(f::MergerForest)     = length(f.trees)
Base.getindex(f::MergerForest, i::Integer) = f.trees[i]
Base.iterate(f::MergerForest, state...)    = iterate(f.trees, state...)
Base.firstindex(f::MergerForest) = firstindex(f.trees)
Base.lastindex(f::MergerForest)  = lastindex(f.trees)

"""
    build_merger_forest(treedir; forest_id=nothing, root_id=nothing,
                        ntrees_max=8, min_halos=1,
                        particle_mass=CT_DEFAULT_PARTICLE_MASS,
                        verbose=true) -> MergerForest

Read the trees of a single consistent-trees **forest** from a rockstar/ctrees
`trees/` directory (the one holding `forests.list`, `locations.dat` and the
`tree_*.dat` files) and return them as a [`MergerForest`](@ref).

Selecting the forest:

* `forest_id`  — read this specific `ForestID`; or
* `root_id`    — read the forest that contains this `TreeRootID`; or
* neither      — **default**: the forest with the largest total halo count
                 (the main cluster's environment).

Only the `ntrees_max` largest trees in the forest (by halo count, a proxy for
mass) with at least `min_halos` halos are read — a forest can contain tens of
thousands of tiny trees, but a merger plot only needs the few massive ones.
Pass `ntrees_max = nothing` to read all of them.

Trees are read by seeking directly to each one via `locations.dat`, so cost
scales with the trees actually requested, not with the full file size.

# Example
```julia
using GadgetIO, MakieHelper, CairoMakie
treedir = "A2256/cr_mhd2/rockstar/trees"
forest  = build_merger_forest(treedir; ntrees_max = 6)
amap    = scale_factor_map(joinpath(treedir, "tree_0_0_0.dat"))
plot_merger_forest(forest;
                   snap_to_scale = s -> amap[s],
                   mass_of       = h -> log10(h.M_TopHat))
```
"""
function build_merger_forest(treedir::AbstractString;
                             forest_id::Union{Integer,Nothing} = nothing,
                             root_id::Union{Integer,Nothing}   = nothing,
                             ntrees_max::Union{Integer,Nothing} = 8,
                             min_halos::Integer = 1,
                             particle_mass::Real = CT_DEFAULT_PARTICLE_MASS,
                             verbose::Bool = true)

    forests_path = joinpath(treedir, "forests.list")
    locs_path    = joinpath(treedir, "locations.dat")
    isfile(forests_path) || error("forests.list not found in $treedir")
    isfile(locs_path)    || error("locations.dat not found in $treedir")

    root2forest = read_forests_list(forests_path)
    locs        = read_tree_locations(locs_path)

    # ----- resolve which forest to read -----
    fid =
        if forest_id !== nothing
            Int64(forest_id)
        elseif root_id !== nothing
            haskey(root2forest, Int64(root_id)) ||
                error("root_id $root_id not present in forests.list")
            root2forest[Int64(root_id)]
        else
            total = Dict{Int64,Int}()       # forest -> total halos
            for (rid, f) in root2forest
                loc = get(locs, rid, nothing)
                loc === nothing && continue
                total[f] = get(total, f, 0) + loc.nhalos
            end
            isempty(total) &&
                error("could not determine a forest (empty forests/locations tables)")
            _, fbig = findmax(total)        # (max total halos, its ForestID)
            fbig
        end

    # ----- roots of this forest that we have locations for, largest first ----
    roots = Int64[rid for (rid, f) in root2forest if f == fid && haskey(locs, rid)]
    isempty(roots) && error("forest $fid has no trees listed in locations.dat")
    filter!(rid -> locs[rid].nhalos >= min_halos, roots)
    sort!(roots; by = rid -> locs[rid].nhalos, rev = true)
    if ntrees_max !== nothing && length(roots) > ntrees_max
        roots = roots[1:ntrees_max]
    end

    verbose && @info @sprintf("forest %d: reading %d of its trees (largest = %d halos)",
                              fid, length(roots), locs[roots[1]].nhalos)

    # ----- read each selected tree by direct seek; open each file once -----
    trees      = Vector{MergerTree}(undef, length(roots))
    headers    = Dict{String,Dict{String,Int}}()
    open_files = Dict{String,IOStream}()
    # read in (file, offset) order for sequential disk access
    order = sortperm(roots; by = rid -> (locs[rid].filename, locs[rid].offset))
    try
        for k in order
            loc   = locs[roots[k]]
            fpath = joinpath(treedir, loc.filename)
            io    = get!(() -> open(fpath, "r"), open_files, loc.filename)
            cols  = get!(() -> _read_ct_header(fpath), headers, loc.filename)
            trees[k] = _read_ct_tree_at(io, loc.offset, cols; particle_mass)
        end
    finally
        for io in values(open_files)
            close(io)
        end
    end

    return MergerForest(fid, trees, roots)
end
