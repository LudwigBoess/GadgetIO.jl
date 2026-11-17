"""
High-level traversal helpers on top of `MergerTree`.

These are written for visualisation / analysis users — they convert the
on-disk 0-based linked-list representation into convenient Julia
iterators and derived per-halo quantities (Num_prog, MMP, Coprog,
log10 mass, ...).
"""

"""
    progenitors(tree, i_julia) -> Vector{Int}

1-based Julia indices of all direct progenitors of halo `i_julia`,
ordered with the most massive (MMP) first — i.e. the order in which the
linked list `FirstProgenitor → NextProgenitor → …` was constructed by
`set_progenitor_pointers!`.
"""
function progenitors(tree::MergerTree, i_julia::Integer)
    result = Int[]
    h = tree[i_julia]
    p = h.FirstProgenitor
    while p >= 0
        push!(result, p + 1)
        p = tree[p + 1].NextProgenitor
    end
    return result
end


"""
    num_progenitors(tree, i_julia) -> Int

Direct progenitor count of halo `i_julia` (1-based).
"""
num_progenitors(tree::MergerTree, i_julia::Integer) =
    length(progenitors(tree, i_julia))


"""
    is_mmp(tree, i_julia) -> Bool

`true` if halo `i_julia` is the most-massive progenitor of its
descendant (i.e. it sits at the head of the descendant's
`FirstProgenitor` chain).
"""
function is_mmp(tree::MergerTree, i_julia::Integer)
    h = tree[i_julia]
    h.Descendant < 0 && return false
    d = tree[h.Descendant + 1]
    return d.FirstProgenitor + 1 == i_julia
end


"""
    coprogenitor_id(tree, i_julia) -> Union{Int,Nothing}

For a halo at the head of a progenitor chain (`is_mmp == true`), return
the next progenitor in the same chain (its "coprogenitor") as a 1-based
index, or `nothing` if there is none. Mirrors galaxybox's `Coprog_ID`
column.
"""
function coprogenitor_id(tree::MergerTree, i_julia::Integer)
    h = tree[i_julia]
    return h.NextProgenitor < 0 ? nothing : Int(h.NextProgenitor + 1)
end


"""
    root_index(tree) -> Int

Index (1-based) of the tree's root halo: the unique halo with
`Descendant == -1`. The depth-first emit order from `generate_trees`
places it at position 1, but we recompute defensively.
"""
function root_index(tree::MergerTree)
    for i in eachindex(tree.halos)
        tree.halos[i].Descendant < 0 && return i
    end
    error("MergerTree has no root (no halo with Descendant == -1)")
end


"""
    walk_main_branch(tree, i_julia=root_index(tree)) -> Vector{Int}

Walk from a halo backwards in time along its most-massive progenitor
branch. Returns 1-based indices ordered from the starting halo
backwards.
"""
function walk_main_branch(tree::MergerTree, i_julia::Integer = root_index(tree))
    chain = Int[i_julia]
    while true
        h = tree[chain[end]]
        h.FirstProgenitor < 0 && break
        push!(chain, h.FirstProgenitor + 1)
    end
    return chain
end


"""
    snap_to_scale_default(snap, last_snap)

Trivial default scale-factor for visualisation when the caller does not
supply an explicit `snap → a` mapping: `a = (snap + 1) / (last_snap + 1)`.
"""
snap_to_scale_default(snap::Integer, last_snap::Integer) =
    (snap + 1) / (last_snap + 1)
