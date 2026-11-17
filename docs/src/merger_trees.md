# Merger Trees

`GadgetIO.jl` ports the original Springel/Lemson [L-BaseTree](https://wwwmpa.mpa-garching.mpg.de/galform/agnpaper/) /
L-HaloTrees C codes (as used for the Millennium and Magneticum projects)
to Julia. The implementation reads/writes the same binary file format as
the C reference, so existing `trees_NNN.X` files can be loaded directly,
and trees built from Julia can be consumed by any L-HaloTrees-compatible
downstream tool.

Two construction stages are involved:

| Stage         | C code        | Julia entry point                |
| :------------ | :------------ | :------------------------------- |
| per-snapshot descendant matching | `L-BaseTree`  | [`build_basetree`](@ref)         |
| stitching into trees             | `L-HaloTrees` | [`build_halotree`](@ref)         |

You'll typically only call [`build_basetree`](@ref) once per consecutive
snapshot pair (to produce a `sub_desc_NNN` file each), then a single
[`build_halotree`](@ref) call that consumes all of them.

## Data structures

The on-disk halo node is a 104-byte isbits struct that matches the
`halo_data` C struct exactly. See the API reference for full field
listings of [`MergerTreeHalo`](@ref), [`MergerTree`](@ref) and
[`MergerTreeFile`](@ref).

Pointer fields (`Descendant`, `FirstProgenitor`, `NextProgenitor`,
`FirstHaloInFOFgroup`, `NextHaloInFOFgroup`) are **0-based indices**
into the tree's halo array (matching the on-disk convention). A value
of `-1` means "no halo". The sentinel is exported as [`NO_HALO`](@ref):

```@docs
NO_HALO
```

When walking trees from Julia, convert to a 1-based index with

```@docs
treeidx
```

## Reading existing trees

L-HaloTrees writes one binary file per "tree file" (typically one per
sub-file of the underlying simulation):

```
Int32  ntrees
Int32  tot_nhalos
Int32  nhalos_per_tree[ntrees]
MergerTreeHalo halos[tot_nhalos]      # depth-first, concatenated
```

To read a single such file:

```@docs
read_merger_tree_file
```

To load every `trees_<snapnum>.X` in a directory at once:

```@docs
read_merger_trees
```

To write back:

```@docs
write_merger_tree_file
```

The companion per-snapshot descendant files produced by L-BaseTree have
their own simple format ([`SubDesc`](@ref)) and are read/written via

```@docs
read_sub_desc
write_sub_desc
```

## Building trees from subfind output

### Step 1 — descendants per snapshot

For two consecutive subfind catalogues `A` (source snap) and `B`
(target snap), [`build_basetree`](@ref) finds the descendant of each
subhalo in `A` by maximising the rank-weighted ID-overlap score

```
score(B_j | A_i) = Σ_{p ∈ A_i ∩ B_j}  1 / (rank_A(p) + 1)^α
```

with `α = 2/3`. Rank weighting biases the match toward the most-bound
particles, which are robust under stripping. An optional `C` catalogue
(snap+2) enables the **decide_upon_descendant** re-routing pass: if
the primary `B`-descendant is part of a merger but a clean `C`-descendant
exists, the halo is redirected to `C`.

```@docs
build_basetree
```

The result is written to disk as `sub_desc_NNN`:

```julia
desc = build_basetree(sub_base_A, sub_base_B, snapnum_B;
                       sub_base_C, snapnum_C = snapnum_B + 1)
write_sub_desc("groups_NNN/sub_desc_NNN", desc)
```

Lower-level building blocks are available if you need finer control
(operate on the [`SubhaloCatalogue`](@ref) container directly):

```@docs
load_subhalo_catalogue
determine_descendants
count_progenitors
decide_upon_descendant!
```

### Step 2 — stitching into trees

Once every `sub_desc_NNN` exists, [`build_halotree`](@ref) loads all
subfind catalogues and descendant files (one [`SubfindMeta`](@ref)
per snapshot), links FOF groups and progenitor chains, then walks each
root subhalo in depth-first order to emit a [`MergerTreeFile`](@ref):

```@docs
load_subfind_meta
assemble_global_halos
set_progenitor_pointers!
generate_trees
build_halotree
```

The end-to-end pipeline reduces to a few calls:

```julia
metas       = Dict(s => load_subfind_meta(sub_base(s)) for s in first_snap:last_snap)
descendants = Dict(s => read_sub_desc(sub_desc_path(s)) for s in first_snap:last_snap-1)

treefile = build_halotree(metas, descendants, first_snap, last_snap)
write_merger_tree_file("treedata/trees_$(last_snap).0", treefile)
```

If your simulation uses non-default overdensity-mass block names,
pass them through `load_subfind_meta`:

```julia
load_subfind_meta(sub_base; mass1 = "MTOP", mass2 = "M5CC", mass3 = "M500")
```

The defaults (`MVIR`, `M5CC`, `M500`) mirror the L-HaloTrees C macros
`MASS1_STRING`, `MASS2_STRING`, `MASS3_STRING`.

## Walking trees

Once you have a [`MergerTree`](@ref), the helper functions below give
1-based Julia indices over the depth-first halo array:

```@docs
root_index
walk_main_branch
progenitors
num_progenitors
is_mmp
coprogenitor_id
```

A typical analysis pattern — pulling the main-branch mass history of
the largest tree in a file:

```julia
treefile  = read_merger_tree_file("treedata/trees_063.0")
tree      = treefile[argmax(length.(treefile.trees))]
chain     = walk_main_branch(tree)                    # 1-based indices
masses    = [tree[i].M_Crit200 for i in chain]
snaps     = [tree[i].SnapNum   for i in chain]
```

## Visualisation

The galaxybox-style merger tree plot
(see [galaxybox `TemporalTreePlotter`](https://github.com/jaoleary/galaxybox/blob/main/src/galaxybox/visualization/tree.py))
is implemented in the companion package
[MakieHelper.jl](https://github.com/LudwigBoess/MakieHelper.jl). It
ships:

* `tree_layout(tree::MergerTree; …)` — pure-Julia port of the
  centred/linear layout algorithm. Returns a `MergerTreeLayout` with
  `(x, y, mass, edges)` ready to hand to any plotting backend.
* `plot_merger_tree(tree; …)` — Makie-native renderer that produces a
  galaxybox-style figure (edges in black, nodes coloured by mass).

```julia
using GadgetIO, MakieHelper, CairoMakie

f    = read_merger_tree_file("treedata/trees_063.0")
tree = f[argmax(length.(f.trees))]

fig, ax = plot_merger_tree(tree;
    mass_of   = h -> log10(h.M_Crit200),
    max_mu    = 1e3,
    plot_style = :centered,
)
save("merger_tree.png", fig)
```

Useful keyword arguments:

| Keyword         | Default                                  | Effect                                                                |
| :-------------- | :--------------------------------------- | :-------------------------------------------------------------------- |
| `snap_to_scale` | `(s+1)/(last_snap+1)`                    | maps `SnapNum` → scale factor for the x-axis                          |
| `mass_of`       | `h -> log10(max(Int(h.Len), 1))`         | per-halo scalar driving the colour                                    |
| `min_scale`     | `0.0`                                    | halos below this scale are dropped from the layout                    |
| `max_mu`        | `10_000.0`                               | mass-ratio cut between MMP and a coprogenitor; tighter ⇒ cleaner plot |
| `plot_style`    | `:centered`                              | `:centered` averages a coprog's y with the MMP's; `:linear` stacks    |
| `colormap`      | `:jet`                                   | Makie colormap symbol                                                 |
| `vmin` / `vmax` | autoscale from valid masses              | colorbar limits                                                       |

For real `snap → a` mapping, build a `Dict{Int,Float64}` from each
snapshot's header:

```julia
a_of_snap = Dict(s => Float64(read_header(snap_base(s)).time)
                 for s in first_snap:last_snap)
plot_merger_tree(tree; snap_to_scale = s -> a_of_snap[s], ...)
```

## Tips and gotchas

* **Subfind block parttype**. `SLEN`, `SOFF` and other subhalo blocks
  live on the subhalo parttype; `PID` lives on the ID parttype.
  [`load_subhalo_catalogue`](@ref) discovers both automatically and
  issues one batched read per parttype.
* **SOFF flavour auto-detection**. Different SUBFIND variants write
  `SOFF` either already-global (`NEW_SUBFIND` without `LGADGET3`) or
  local-to-each-sub-file (`LGADGET3`-style). [`load_subhalo_catalogue`](@ref)
  picks the interpretation that produces a consistent
  `(suboffset + sublen) ≤ length(ids)` layout, so callers need not know
  which compile flags built their subfind output.
* **Threading**. [`determine_descendants`](@ref) and
  [`read_subfind_blocks`](@ref) both pick up extra cores automatically.
  Launch Julia with `julia -t N` to scale; per-task scratch buffers
  are allocated inside [`determine_descendants`](@ref) so the result is
  deterministic and safe under Julia ≥ 1.10 task migration.
* **Reading vs. building**. If a downstream tool (e.g. a SAM) wrote
  L-HaloTrees files already, you can skip [`build_basetree`](@ref) /
  [`build_halotree`](@ref) entirely and start from
  [`read_merger_tree_file`](@ref).
* **Pointer indexing**. On-disk and in-memory pointer fields are
  0-based with `-1` as the "no halo" sentinel, matching the C/IDL
  convention. Use [`treeidx`](@ref) to convert to a 1-based Julia
  index when indexing into `tree.halos`.

## Batched subfind reads (advanced)

The merger-tree pipeline relies on a single helper that reads many
subfind blocks in one pass over the sub-files, caching headers, info
and block positions, and threading the IO. You can use it directly if
you need to pull many blocks at once for any other reason:

```@docs
read_subfind_blocks
```
