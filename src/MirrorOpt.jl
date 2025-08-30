using JuMP
using HiGHS  # or GLPK, Cbc, Gurobi, CPLEX

if !isdefined(Main, :CONSTANTS_INCLUDED)
    include("Constants.jl")
    const CONSTANTS_INCLUDED = true
end

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
    # seen keys are tuples of variable length (one entry per cell)
    seen = Set{Tuple{Vararg{Tuple{Int,Int}}}}()
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

# Each edge is identified by (cell, dir), with dir ∈ (:U,:D,:L,:R).
# Compute the set of external edges of a placed shape (grid aligned).
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
        # For :D the mirror line is y = ey + 0.5 (neighbor below), for :U it's y = ey - 0.5
        δ = (d == :D) ? 1 : -1
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
            # Only allow original anchor points (even coordinates) to avoid starting on padded cells
            for ax in 1:2:(W - maxx)
                for ay in 1:2:(H - maxy)
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
    # seen keys: pair of variable-length tuples for cu and cl (supports expanded shapes)
    seen = Set{Tuple{Tuple{Vararg{Tuple{Int,Int}}}, Tuple{Vararg{Tuple{Int,Int}}}}}()
        for U in ups[t]
            for e in external_edges(U)
                L = reflect_across_edge(U, e)
                if all(in_bounds, L)
                    cu = sort(U)
                    cl = sort(L)
                    key = (Tuple(cu), Tuple(cl))
                    if !(key in seen)
                        push!(seen, key)
                        # illegal if the mirrored shape shares any cells with the original
                        if !isempty(intersect(Set(cu), Set(cl)))
                            continue
                        end
                        wU = sum(A[y, x] for (x,y) in cu)
                        wL = sum(A[y, x] for (x,y) in cl)
                        total_sum = wU - wL
                        # Just a heuristic, risks unoptimality, skipping negative placements
                        if total_sum >= 0 # TODO: try without 
                            push!(plist, (cellsU = cu, cellsL = cl, weight = total_sum))
                        end
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
    # Flatten to indexable list
    index = Dict{Symbol, Vector{Int}}()
    P = Vector{NamedTuple}()  # each element has fields: t, cells, weight
    for t in LETTERS
        index[t] = Int[]
        for pr in pairs[t]
            push!(index[t], length(P)+1)
            cellsAll = union(Set(pr.cellsU), Set(pr.cellsL)) |> collect
            # store both the combined cells and the original U/L partitions
            push!(P, (t=t, cells=cellsAll, cellsU = pr.cellsU, cellsL = pr.cellsL, weight=pr.weight))
        end
    end

    # Build an W x H matrix where each cell contains the list of pair indices
    # that would occupy that cell or any of its 8-neighbors (used for cell-wise <=1 constraints)
    occ8 = [Int[] for x in 1:W, y in 1:H]

    # helper for 8-neighborhood including the cell itself
    function neighborhood8((x,y))
        ((x,y),
         (x-1,y-1),(x,y-1),(x+1,y-1),
         (x-1,y),(x+1,y),
         (x-1,y+1),(x,y+1),(x+1,y+1))
    end

    for (pid, pr) in enumerate(P)
        # For each cell touched by the pair, add pid to that cell and its neighbors
        for c in pr.cells
            for nb in neighborhood8(c)
                if in_bounds(nb)
                    nx, ny = nb
                    push!(occ8[nx, ny], pid)
                end
            end
        end
    end

    return P, index, occ8
end

# Build and solve the MILP
function solve_pentomino(A)
    ups   = enumerate_uppercase_placements()
    pairs = build_pairs(A, ups)
    println("DEBUG: built $(sum(length(pairs[t]) for t in LETTERS)) pairs")
    P, index, occ8 = build_conflicts(pairs)

    # If no pairs were generated, return early (no model to build)
    if length(P) == 0
        println("WARNING: no candidate pairs; skipping model build")
        return :NO_PAIRS, NaN, P, Int[]
    end

    # Build and solve model inside try/catch to avoid throwing when solver/setup fails
    status = :ERROR
    obj = NaN
    chosen = Int[]
    try
        model = Model(HiGHS.Optimizer)  # or GLPK.Optimizer / Cbc.Optimizer / Gurobi.Optimizer
        #set_optimizer_attribute(model, "time_limit", 1800.0)  # 30-minute time limit
        #set_silent(model)

        @variable(model, y[1:length(P)], Bin)

    # at most one pair per letter (allow skipping letters for speed/score)
    for t in LETTERS
        ids = index[t]
        if !isempty(ids)
            #@constraint(model, sum(y[i] for i in ids) <= 1)
            @constraint(model, sum(y[i] for i in ids) == 1)
        end
    end
    # no overlap is already implied by the conflicts set built from overlap;
    # but adding cell-wise ≤1 can strengthen. Build once:
    # (Optional) Build cell-wise cap tightening:
    # ...

    # Enforce non-overlap / non-touching by cell-wise constraints using occ8
    # Deduplicate indices per cell, avoid shadowing JuMP variable `y` by using xi, yi.
    for xi in 1:W, yi in 1:H
        ids = unique(occ8[xi, yi])
        if length(ids) > 1
            @constraint(model, sum(y[i] for i in ids) ≤ 1)
        end
    end

        @objective(model, Max, sum(P[i].weight * y[i] for i in 1:length(P)))

        optimize!(model)

        status = termination_status(model)
        # Safely attempt to read objective and variable values; handle cases with no solution
        try
            obj    = objective_value(model)
            chosen = findall(i -> value(y[i]) > 0.5, 1:length(P))
        catch _err
            obj = NaN
            chosen = Int[]
        end
    catch err
        # Log and return a safe failure result instead of throwing
        println("ERROR: building or solving model failed: ", err)
        return :ERROR, NaN, P, Int[]
    end
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
            # mark uppercase and mirrored (lowercase) cells distinctly
            up = first(String(P[i].t))
            low = lowercase(string(up))[1]
            for (x,y) in P[i].cellsU
                grid[y,x] = up
            end
            for (x,y) in P[i].cellsL
                grid[y,x] = low
            end
        end
        # Print a reduced grid sampling the top-left cell of each 2x2 block
        reduced_h = div(H, 2)
        reduced_w = div(W, 2)
        println("solution grid (reduced):")
        for ry in 1:reduced_h
            rowchars = Vector{Char}(undef, reduced_w)
            for rx in 1:reduced_w
                srcx = 2*(rx-1) + 1
                srcy = 2*(ry-1) + 1
                rowchars[rx] = grid[srcy, srcx]
            end
            println(join(rowchars))
        end
    end
end

return 0