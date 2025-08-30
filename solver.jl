using JuMP
using HiGHS  # or GLPK, Cbc, Gurobi

const N = 16

const POINT_VALUE_MATRIX = [
    4 2 6 3 1 5 8 7 6 9 1 9 4 7 8 3;
    8 7 1 5 6 3 2 4 9 5 7 3 6 9 1 8;
    8 9 5 8 4 7 9 6 1 8 3 4 8 6 5 7;
    1 9 3 2 7 6 1 2 4 9 5 7 3 2 8 6;
    9 7 4 6 5 1 7 8 6 3 2 8 1 5 2 9;
    2 8 4 1 9 9 3 7 2 6 4 5 8 3 9 1;
    3 5 8 5 2 6 4 8 9 1 6 1 2 4 3 9;
    5 1 8 4 3 9 2 1 7 5 3 7 9 7 4 2;
    7 8 2 9 6 4 7 3 5 4 8 9 6 1 7 5;
    4 2 6 3 1 5 9 5 8 2 9 4 3 6 1 8;
    5 4 1 6 9 2 3 7 1 5 4 6 7 2 8 3;
    6 3 7 2 4 8 5 6 9 3 2 7 1 8 5 4;
    2 9 5 1 7 4 1 3 6 5 8 4 7 3 2 6;
    7 4 3 5 2 6 6 1 2 4 1 7 5 9 3 8;
    1 6 2 4 8 5 3 7 4 9 3 2 6 1 7 5;
    7 5 9 3 6 9 8 2 7 4 6 2 8 5 8 1
]
# TODO: verify value matrix matches official problem statement

const EXAMPLE_SOLUTION = [
    '.' 'f' '.' 'v' 'v' 'v' 'V' 'V' 'V' '.' '.' '.' 'z' 'Z' '.' '.';
    'f' 'f' '.' '.' '.' 'v' 'V' '.' '.' '.' 'z' 'z' 'z' 'Z' 'Z' 'Z';
    '.' 'f' 'f' '.' '.' 'v' 'V' '.' '.' '.' 'z' '.' '.' '.' '.' 'Z';
    '.' 'F' 'F' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.';
    'F' 'F' '.' '.' 't' '.' '.' 'T' '.' '.' '.' 'X' '.' '.' 'x' '.';
    '.' 'F' '.' '.' 't' '.' '.' 'T' '.' '.' 'X' 'X' 'X' 'x' 'x' 'x';
    '.' '.' '.' 't' 't' 't' 'T' 'T' 'T' '.' '.' 'X' '.' '.' 'x' '.';
    'I' 'i' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.' '.';
    'I' 'i' '.' 'L' '.' '.' '.' '.' '.' '.' 'N' 'N' '.' 'U' 'U' 'U';
    'I' 'i' '.' 'L' 'L' 'L' 'L' '.' 'N' 'N' 'N' '.' '.' 'U' '.' 'U';
    'I' 'i' '.' 'l' 'l' 'l' 'l' '.' 'n' 'n' 'n' '.' '.' 'u' '.' 'u';
    'I' 'i' '.' 'l' '.' '.' '.' '.' '.' '.' 'n' 'n' '.' 'u' 'u' 'u';
    '.' '.' '.' '.' '.' 'y' 'y' 'y' 'y' '.' '.' '.' '.' '.' '.' '.';
    'P' 'P' 'p' 'p' '.' '.' '.' 'y' '.' '.' 'w' '.' '.' '.' '.' 'W';
    'P' 'P' 'p' 'p' '.' '.' '.' 'Y' '.' '.' 'w' 'w' '.' '.' 'W' 'W';
    'P' '.' '.' 'p' '.' 'Y' 'Y' 'Y' 'Y' '.' '.' 'w' 'w' 'W' 'W' '.'
]
# TODO: verify, seems fine

const A = POINT_VALUE_MATRIX
const H, W = size(A)
const LETTERS = [:F, :I, :L, :N, :P, :T, :U, :V, :W, :X, :Y, :Z]

# Derive pentomino base shapes from the uppercase placements in EXAMPLE_SOLUTION.
# Shapes are returned as 5 relative (x,y) offsets with min x/y at 0.
const SHAPES = let
    function extract_shape(letter::Char, grid::AbstractMatrix{Char})
        coords = Tuple{Int,Int}[]
        for y in 1:size(grid,1), x in 1:size(grid,2)
            if grid[y,x] == letter
                push!(coords, (x,y))
            end
        end
        length(coords) == 5 || error("Expected 5 cells for '$letter' in EXAMPLE_SOLUTION, got $(length(coords))")
        minx = minimum(first.(coords))
        miny = minimum(last.(coords))
        sort([(x - minx, y - miny) for (x,y) in coords])
    end

    shapes = Dict{Symbol, Vector{Tuple{Int,Int}}}()
    for s in LETTERS
        ch = first(String(s))  # uppercase Char for the symbol
        shapes[s] = extract_shape(ch, EXAMPLE_SOLUTION)
    end
    shapes
end

print(SHAPES)

exit(0)

# 2D transforms
rot90((x,y)) = (-y, x)
rot(p::Tuple{Int,Int}, k::Int) = k == 0 ? p : rot(rot90(p), k-1)
reflX((x,y)) = (-x, y)   # mirror horizontally
reflY((x,y)) = (x, -y)   # mirror vertically

# Normalize a shape: shift so min x and min y are zero; sort for canonical
function normalize(cells::Vector{Tuple{Int,Int}})
    xs = [x for (x,_) in cells]
    ys = [y for (_,y) in cells]
    dx, dy = minimum(xs), minimum(ys)
    sort!([(x-dx, y-dy) for (x,y) in cells])
end

# Generate all unique rotations of a base shape (uppercase side, no mirror)
function generate_rotations(base)
    seen = Set{NTuple{5,Tuple{Int,Int}}}()
    outs = Vector{Vector{Tuple{Int,Int}}}()
    for k in 0:3
        cells = [rot(c,k) for c in base]
        norm = normalize(cells)
        key = Tuple(norm)
        if !(key in seen)
            push!(outs, norm)
            push!(seen, key)
        end
    end
    outs
end

# Translate to absolute grid cells with anchor shift (ax, ay)
translate(cells, ax, ay) = [(x+ax, y+ay) for (x,y) in cells]

in_bounds((x,y)) = 1 ≤ x ≤ W && 1 ≤ y ≤ H

# Compute the set of external edges of a placed shape (grid aligned).
# Each edge is identified by (cell, dir), with dir ∈ (:U,:D,:L,:R).
const DIRS = Dict(:U => (0,1), :D => (0,-1), :L => (-1,0), :R => (1,0))
function external_edges(cells::Vector{Tuple{Int,Int}})
    cellset = Set(cells)
    edges = Vector{Tuple{Tuple{Int,Int},Symbol}}()
    for (x,y) in cells
        for (d, (dx,dy)) in DIRS
            nb = (x+dx, y+dy)
            if !(nb in cellset)
                push!(edges, ((x,y), d))
            end
        end
    end
    edges
end

# Reflect a placed shape across the grid line that coincides with a specific edge.
# For an edge at cell (x,y) with :R, the mirror line is x = x + 0.5 (vertical).
function reflect_across_edge(cells::Vector{Tuple{Int,Int}}, edge::Tuple{Tuple{Int,Int},Symbol})
    ((ex,ey), d) = edge
    if d == :L || d == :R
        # vertical mirror line: x = ex ± 0.5
        δ = (d == :R) ? 1 : -1
        return [(2*ex + δ - x, y) for (x,y) in cells]
    else
        # horizontal mirror line: y = ey ± 0.5
        δ = (d == :U) ? 1 : -1
        return [(x, 2*ey + δ - y) for (x,y) in cells]
    end
end

# Enumerate all uppercase placements (shape rotations and translations)
function enumerate_uppercase_placements()
    ups = Dict{Symbol, Vector{Vector{Tuple{Int,Int}}}}()
    for t in LETTERS
        base = SHAPES[t]
        rots = generate_rotations(base)
        placements = Vector{Vector{Tuple{Int,Int}}}()
        for rshape in rots
            # bounding to keep in the grid (rshape is normalized to min x/y = 0)
            xs = (x for (x,_) in rshape); ys = (y for (_,y) in rshape)
            maxx, maxy = maximum(xs), maximum(ys)
            for ax in 1:(W - maxx)
                for ay in 1:(H - maxy)
                    placed = translate(rshape, ax, ay)  # 1-based grid coords
                    # in-bounds check is redundant given bounds, but keep it safe
                    if all(in_bounds, placed)
                        push!(placements, placed)
                    end
                end
            end
        end
        ups[t] = placements
    end
    ups
end

# Build all valid mirrored-and-touching pairs for each letter
function build_pairs(A, ups)
    # Return:
    # pairs_by_letter::Dict{Symbol, Vector{NamedTuple}} with fields:
    #   cellsU::Vector{Tuple{Int,Int}}, cellsL::Vector{Tuple{Int,Int}}, weight::Int
    pairs = Dict{Symbol, Vector{NamedTuple}}()
    for t in LETTERS
        plist = NamedTuple[]
        seen = Set{Tuple{NTuple{5,Tuple{Int,Int}},NTuple{5,Tuple{Int,Int}}}}()
        for U in ups[t]
            for e in external_edges(U)
                L = reflect_across_edge(U, e)
                if all(in_bounds, L)
                    # ensure L is a pure mirror of U (already the case by construction)
                    # de-dup: canonical order
                    cu = sort(U)
                    cl = sort(L)
                    key = (Tuple(cu), Tuple(cl))
                    if !(key in seen)
                        push!(seen, key)
                        wU = sum(A[y, x] for (x,y) in cu)
                        wL = sum(A[y, x] for (x,y) in cl)
                        push!(plist, (cellsU = cu, cellsL = cl, weight = wU - wL))
                    end
                end
            end
        end
        pairs[t] = plist
    end
    pairs
end

# Build conflicts between pairs of different letters (overlap or 8-neighbor touch)
function build_conflicts(pairs)
    # flatten to indexable list
    index = Dict{Symbol, Vector{Int}}()
    P = Vector{NamedTuple}()  # each element has fields: t, cells, weight
    for t in LETTERS
        index[t] = Int[]
        for pr in pairs[t]
            # Heuristic: drop negative-score placements entirely (y == 0)
            pr.weight < 0 && continue
            push!(index[t], length(P)+1)
            cellsAll = union(Set(pr.cellsU), Set(pr.cellsL)) |> collect
            push!(P, (t=t, cells=cellsAll, weight=pr.weight))
        end
    end
    # cell->indices map for fast overlap
    bycell = Dict{Tuple{Int,Int}, Vector{Int}}()
    for (pid, pr) in enumerate(P)
        for c in pr.cells
            push!(get!(bycell, c, Int[]), pid)
        end
    end
    # helper for 8-neighborhood
    function neighbors8((x,y))
        ((x-1,y-1),(x,y-1),(x+1,y-1),
         (x-1,y),(x+1,y),
         (x-1,y+1),(x,y+1),(x+1,y+1))
    end
    conflicts = Set{Tuple{Int,Int}}()
    # overlaps
    for (_, lst) in bycell
        for i in 1:length(lst), j in (i+1):length(lst)
            p, q = lst[i], lst[j]
            if P[p].t != P[q].t
                push!(conflicts, (min(p,q), max(p,q)))
            end
        end
    end
    # 8-neighbor touches
    occ = Dict{Tuple{Int,Int}, Vector{Int}}()
    for (pid, pr) in enumerate(P)
        for c in pr.cells
            push!(get!(occ, c, Int[]), pid)
        end
    end
    for (pid, pr) in enumerate(P)
        for c in pr.cells
            for nb in neighbors8(c)
                for q in get(occ, nb, Int[])
                    if pid < q && P[pid].t != P[q].t
                        push!(conflicts, (pid, q))
                    end
                end
            end
        end
    end
    return P, index, collect(conflicts)
end

# Build and solve the MILP
function solve_pentomino(A)
    ups   = enumerate_uppercase_placements()
    pairs = build_pairs(A, ups)
    P, index, conflicts = build_conflicts(pairs)

    model = Model(HiGHS.Optimizer)  # or GLPK.Optimizer / Cbc.Optimizer / Gurobi.Optimizer
    set_optimizer_attribute(model, "time_limit", 120.0)  # 2-minute time limit
    set_silent(model)

    @variable(model, y[1:length(P)], Bin)
    # at most one pair per letter (allow skipping letters for speed/score)
    for t in LETTERS
        ids = index[t]
        if !isempty(ids)
            @constraint(model, sum(y[i] for i in ids) <= 1)
        end
    end
    # no overlap is already implied by the conflicts set built from overlap;
    # but adding cell-wise ≤1 can strengthen. Build once:
    # (Optional) Build cell-wise cap tightening:
    # ...

    # no-touch-even-diagonal across different letters
    for (p,q) in conflicts
        @constraint(model, y[p] + y[q] ≤ 1)
    end

    @objective(model, Max, sum(P[i].weight * y[i] for i in 1:length(P)))

    println("WARNING: Skipping optimization!")
    # optimize!(model)

    status = termination_status(model)
    obj    = objective_value(model)
    chosen = findall(i -> value(y[i]) > 0.5, 1:length(P))
    return status, obj, P, chosen
end

# Simple runner when executing this file

function run()
    status, best, P, chosen = solve_pentomino(A)
    println("status: ", status)
    println("objective: ", best)
    println("selected pairs: ", length(chosen))
    if !isempty(chosen)
        grid = fill('.', H, W)
        for i in chosen
            # mark uppercase and lowercase cells distinctly
            up = first(String(P[i].t))
            low = lowercase(string(up))[1]
            # Attempt to recover U/L from cells by symmetry: not stored in P, mark all as up for now
            for (x,y) in P[i].cells
                grid[y,x] = up
            end
        end
        println("solution grid:")
        for y in 1:H
            println(join(grid[y,:]))
        end
    end
end