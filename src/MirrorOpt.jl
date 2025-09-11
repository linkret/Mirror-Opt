using JuMP
using HiGHS  # or GLPK, Cbc, Gurobi, CPLEX
using CPLEX
using Random

using Printf: @sprintf
import MathOptInterface as MOI

# Seed RNG from current time (similar to C++ time(0))
Random.seed!(UInt64(Base.time_ns()))

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

# Generate all unique rotations and reflections of a base shape
function generate_rotations(base)
    # seen keys are tuples of variable length (one entry per cell)
    seen = Set{Tuple{Vararg{Tuple{Int,Int}}}}()
    outs = Vector{Vector{Tuple{Int,Int}}}()
    # iterate over horizontal and vertical mirror options; normalize dedups
    for mx in (false, true)
        for my in (false, true)
            shape0 = base
            if mx
                shape0 = [reflX(p) for p in shape0]
            end
            if my
                shape0 = [reflY(p) for p in shape0]
            end
            for k in 0:3
                cells = [rot(p, k) for p in shape0]
                norm = normalize(cells)
                key = Tuple(norm)
                if !(key in seen)
                    push!(outs, norm)
                    push!(seen, key)
                end
            end
        end
    end
    return outs
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

    # Collect all boundary edges per direction, then keep only extremal ones:
    # L: leftmost x (tie: smallest y)
    # R: rightmost x (tie: smallest y)
    # U: topmost (smallest y; tie: smallest x)
    # D: bottommost (largest y; tie: smallest x)
    cands = Dict{Symbol, Vector{Tuple{Tuple{Int,Int},Symbol}}}(
        :L => Tuple{Tuple{Int,Int},Symbol}[],
        :R => Tuple{Tuple{Int,Int},Symbol}[],
        :U => Tuple{Tuple{Int,Int},Symbol}[],
        :D => Tuple{Tuple{Int,Int},Symbol}[],
    )
    for (x,y) in cells
        for (d, (dx,dy)) in DIRS
            nb = (x+dx, y+dy)
            if !(nb in cellset)
                push!(cands[d], ((x,y), d))
            end
        end
    end

    edges = Tuple{Tuple{Int,Int},Symbol}[]

    # Select one candidate per direction based on extremal rule
    if !isempty(cands[:L])
        best = nothing; bx = typemax(Int); by = typemax(Int)
        for e in cands[:L]
            (x,y) = e[1]
            if x < bx || (x == bx && y < by)
                best = e; bx = x; by = y
            end
        end
        push!(edges, best)
    end
    if !isempty(cands[:R])
        best = nothing; bx = -typemax(Int); by = typemax(Int)
        for e in cands[:R]
            (x,y) = e[1]
            if x > bx || (x == bx && y < by)
                best = e; bx = x; by = y
            end
        end
        push!(edges, best)
    end
    if !isempty(cands[:U])
        best = nothing; by = typemax(Int); bx = typemax(Int)
        for e in cands[:U]
            (x,y) = e[1]
            if y < by || (y == by && x < bx)
                best = e; by = y; bx = x
            end
        end
        push!(edges, best)
    end
    if !isempty(cands[:D])
        best = nothing; by = -typemax(Int); bx = typemax(Int)
        for e in cands[:D]
            (x,y) = e[1]
            if y > by || (y == by && x < bx)
                best = e; by = y; bx = x
            end
        end
        push!(edges, best)
    end

    return edges
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

# helper for 8-neighborhood excluding the cell itself
function neighbors8((x,y))
    ((x-1,y-1),(x,y-1),(x+1,y-1),(x-1,y),(x+1,y),(x-1,y+1),(x,y+1),(x+1,y+1))
end

function inflate(pr)
    s = union(pr.cellsU, pr.cellsL)
    # add neighbors
    for (x,y) in collect(s)
        for nb in neighbors8((x,y))
            if in_bounds(nb)
                push!(s, nb)
            end
        end
    end
    return s
end

# Enumerate all uppercase placements (shape rotations and translations)
function enumerate_uppercase_placements()
    ups = Dict{Symbol, Vector{Vector{Tuple{Int,Int}}}}()
    for t in LETTERS
        base = SHAPES[t]
        rots = generate_rotations(base)
        placements = Vector{Vector{Tuple{Int,Int}}}()
        for rshape in rots
            # Only allow original anchor points (even coordinates) to avoid starting on padded cells
            for ax in 1-N:2:N+N
                for ay in 1-N:2:N+N
                    placed = translate(rshape, ax, ay)  # 1-based grid coords
                    
                    if !all(in_bounds, placed)
                        continue
                    end
                    
                    # total_sum = sum(A[y, x] for (x,y) in placed)

                    push!(placements, placed)
                end
            end
        end
        ups[t] = placements
    end
    ups
end

# Pretty-print a reduced grid (top-left of each 2x2 block) for selected pairs
function print_solution_grid(P::Vector{NamedTuple}, selected::Vector{Int}; io::IO=stdout)
    if isempty(selected)
        println(io, "(no pairs)")
        return
    end
    grid = fill('.', H, W)
    for i in selected
        up = first(String(P[i].t))
        low = lowercase(string(up))[1]
        for (x,y) in P[i].cellsU
            grid[y,x] = up
        end
        for (x,y) in P[i].cellsL
            grid[y,x] = low
        end
    end
    reduced_h = div(H, 2)
    reduced_w = div(W, 2)
    for ry in 1:reduced_h
        rowchars = Vector{Char}(undef, reduced_w)
        for rx in 1:reduced_w
            srcx = 2*(rx-1) + 1
            srcy = 2*(ry-1) + 1
            rowchars[rx] = grid[srcy, srcx]
        end
        println(io, join(rowchars))
    end
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
                        # Check against example solution to mark warm-start pair
                        upch = first(String(t))
                        loch = lowercase(string(upch))[1]
                        matches_example = all(EXAMPLE_SOLUTION[y, x] == upch for (x,y) in cu) &&
                                          all(EXAMPLE_SOLUTION[y, x] == loch for (x,y) in cl)
                        minx = minimum(x for (x,y) in cu)
                        miny = minimum(y for (x,y) in cu)
                        maxx = maximum(x for (x,y) in cu)
                        maxy = maximum(y for (x,y) in cu)
                        # Heuristic speed-up: skip negative placements
                        if matches_example || total_sum >= 5
                            push!(plist, (
                                cellsU = cu,
                                cellsL = cl,
                                weight = total_sum,
                                is_ws = matches_example,
                                bbox = (minx, miny, maxx, maxy)
                            ))
                        end
                    end
                end
            end
        end

        function halve(cells)
            s = Set{Tuple{Int,Int}}()
            for (x,y) in cells
                push!(s, (x÷2,y÷2))
            end
            return s
        end

        plist = sort(plist, by = x -> (x.bbox, -x.weight))
        inflated = [halve(inflate(p)) for p in plist]
        # Remove duplicates by inflated shape, keeping only one per bbox, largest Weighted

        keep = trues(length(plist))
        
        for ip in eachindex(plist)
            if !keep[ip]
                continue
            end
            p = plist[ip]
            for ir in (ip+1):length(plist)
                r = plist[ir]
                
                if !keep[ir]
                    continue
                end
                
                if p.bbox != r.bbox
                    break
                end

                if inflated[ip] == inflated[ir]
                    keep[ir] = false
                end
            end
        end

        plist = plist[keep]
        plist = sort(plist, by = x -> -x.weight)
        pairs[t] = plist
    end
    pairs
end

# Build conflicts between pairs of different letters (overlap or 8-neighbor touch)
function build_conflicts(pairs)
    # Flatten to indexable list
    index = Dict{Symbol, Vector{Int}}()
    P = Vector{NamedTuple}()  # each element has fields: t, cells, cellsU, cellsL, weight
    WS = Int[]  # flattened y indices that match the example solution
    
    for t in LETTERS
        index[t] = Int[]
        for pr in pairs[t]
            push!(index[t], length(P)+1)
            cellsAll = union(Set(pr.cellsU), Set(pr.cellsL)) |> collect
            # store both the combined cells and the original U/L partitions
            push!(P, (t=t, cells=cellsAll, cellsU = pr.cellsU, cellsL = pr.cellsL, weight=pr.weight))
            if hasproperty(pr, :is_ws) && pr.is_ws
                push!(WS, length(P))
            end
        end
    end

    # Build an W x H matrix where each cell contains the list of pair indices
    # that occupy that cell (NO neighbor expansion here)
    occ = [Int[] for x in 1:W, y in 1:H]

    # populate occ with pair indices for the cells they actually occupy
    for (pid, pr) in enumerate(P)
        for c in pr.cells
            if in_bounds(c)
                nx, ny = c
                push!(occ[nx, ny], pid)
            end
        end
    end

    # build conflicts set: overlaps (same cell) and neighbor touches
    # conflicts = Set{Tuple{Int,Int}}()
    ## overlaps: any two different-letter pairs sharing a cell
    # for x in 1:W, y in 1:H
    #     lst = occ[x,y]
    #     for i in 1:length(lst), j in (i+1):length(lst)
    #         p, q = lst[i], lst[j]
    #         if P[p].t != P[q].t
    #             push!(conflicts, (min(p,q), max(p,q)))
    #         end
    #     end
    # end

    # println("DEBUG: checking neighbor touches")

    # prepare per-pair neighbor-cell lists (this is what the caller expects)
    neighcells = [Vector{Tuple{Int,Int}}() for _ in 1:length(P)]

    for (pid, pr) in enumerate(P)
        # compute the outer-edge neighbor cells of the whole shape (unique)
        neighs = Set{Tuple{Int,Int}}()
        for c in union(pr.cellsU, pr.cellsL)
            for nb in neighbors8(c)
                if in_bounds(nb)
                    push!(neighs, nb)
                end
            end
        end
        # remove cells that are actually occupied by the pair to get the "outer edge"
        for c in union(pr.cellsU, pr.cellsL)
            if c in neighs
                delete!(neighs, c)
            end
        end

        # record neighbor cells for this pair (caller uses this)
        neighcells[pid] = collect(neighs)

        # # iterate outer neighbors (unique) once and add conflicts with pairs occupying them
        # for nb in neighs
        #     nx, ny = nb
        #     for q in occ[nx, ny]
        #         if pid < q && P[pid].t != P[q].t
        #             push!(conflicts, (pid, q))
        #         end
        #     end
        # end
    end

    return P, index, occ, neighcells, WS
end

# Greedy pre-selection (in-place): pick k highest-weight pairs of distinct letters without overlaps.
# Assumes pairs[t] is sorted descending by weight. Mutates `pairs` to prune incompatible pairs and
# keep at most one (chosen) pair per selected letter. Returns the flattened indices (over the mutated order)
# of the chosen pairs.
function pick_greedy!(k::Int, pairs::Dict{Symbol, Vector{NamedTuple}}, randomize::Bool)
    chosen_keys = Vector{Tuple{Symbol, Vector{Tuple{Int,Int}}, Vector{Tuple{Int,Int}}}}()
    chosen_letters = Set{Symbol}()
    # blocked contains all cells occupied by chosen pairs AND their 8-neighborhoods
    blocked = Set{Tuple{Int,Int}}()
    pos = Dict{Symbol, Int}(t => 1 for t in LETTERS)  # retained (not used for scoring now)

    # Weighted "cells taken" cost: occupied cells count as 1; empty neighbor cells count as 1,
    # but if a neighbor is on the outer edge band (first/last two rows/cols), it counts as 2.
    # outer_edge((x,y)) = (x == 1 || x == 2 || x == W-1 || x == W || y == 1 || y == 2 || y == H-1 || y == H)
    outer_edge((x,y)) = (min(x, y) <= 2 || max(x, y) >= N-1)

    function taken_cost(pr)::Int
        # return 1 # TODO: remove
        occ = union(pr.cellsU, pr.cellsL)
        # neighbors excluding occupied
        nbs = setdiff(inflate(pr), occ)
        cost = length(occ)
        for nb in nbs
            cost += outer_edge(nb) ? 2 : 1
        end
        return cost
    end

    overlaps_blocked = function (pr)
        for c in inflate(pr)
            if c in blocked
                return true
            end
        end
        return false
    end

    is_edge_ok = function (t::Symbol, pr)
        # For letter I, require at least one cell on the matrix edge; else reject for greedy
        return true
        if t in [:W]
            return true # TODO: remove, hack
        end
        for (x,y) in union(pr.cellsU, pr.cellsL)
            if x == 1 || x == W || y == 1 || y == H
                #if t in [:P]
                #    return false
                #end
                #if t in [:P]
                #    return false
                #end
                return true
            end
        end
        #if t in [:P]
        #    return true
        #end
        return false
    end

    function best_feasible_for(t::Symbol)
        v = get(pairs, t, NamedTuple[])
        best = nothing
        best_metric = -Inf
        best_weight = -typemax(Int)
        for pr in v
            if !overlaps_blocked(pr) && is_edge_ok(t, pr)
                cells_taken = taken_cost(pr)
                if cells_taken == 0
                    continue
                end
                m = pr.weight / cells_taken
                if randomize
                    m *= (0.7 + 0.6 * rand())
                end
                if m > best_metric || (m == best_metric && pr.weight > best_weight)
                    best_metric = m
                    best_weight = pr.weight
                    best = pr
                end
            end
        end
        return best
    end

    for _ in 1:k
        best_t = nothing
        best_pr = nothing
        best_w = -Inf
        for t in LETTERS
            if t in chosen_letters
                continue
            end
            pr = best_feasible_for(t)
            if pr === nothing
                continue
            end
            m = pr.weight / taken_cost(pr)
            # if t in [:W, :U, :I]
            #     m *= 2 # TODO: remove, temporary hack
            # end
            if randomize
                m *= (0.7 + 0.6 * rand()) # randomization
            end
            if m > best_w
                best_w = m
                best_t = t
                best_pr = pr
            end
        end
        if best_t === nothing
            break
        end

        push!(chosen_keys, (best_t, best_pr.cellsU, best_pr.cellsL))
        push!(chosen_letters, best_t)
    for c in inflate(best_pr); push!(blocked, c); end

        # Mutate: keep only the chosen for best_t
        pairs[best_t] = [best_pr]
        # For other letters (excluding all already chosen), drop overlapping pairs
    for t in LETTERS
            if t in chosen_letters
                continue
            end
            v = get(pairs, t, NamedTuple[])
            if isempty(v)
                continue
            end
            keep = NamedTuple[]
            for pr in v
        if !overlaps_blocked(pr)
                    push!(keep, pr)
                end
            end
            pairs[t] = keep
            pos[t] = min(get(pos, t, 1), length(keep))
            if pos[t] == 0
                pos[t] = 1
            end
        end
    end

    # Compute flattened indices in the mutated (current) pairs order
    chosen_indices = Int[]
    idx = 0
    for t in LETTERS
        for pr in get(pairs, t, NamedTuple[])
            idx += 1
            for (ct, cu, cl) in chosen_keys
                if t == ct && pr.cellsU == cu && pr.cellsL == cl
                    push!(chosen_indices, idx)
                    break
                end
            end
        end
    end
    return chosen_indices
end

const UPS = enumerate_uppercase_placements()
const PAIRS = build_pairs(A, UPS)

# Build and solve the MILP
function solve_pentomino(A, greedy_k::Int = 0, optimize::Bool = true, randomize::Bool = false)
    # ups = deepcopy(UPS)
    pairs = deepcopy(PAIRS)
    picked = pick_greedy!(greedy_k, pairs, randomize)

    println("DEBUG: built $(sum(length(pairs[t]) for t in LETTERS)) pairs")
    P, index, occ, neighcells, WS = build_conflicts(pairs)
    println("DEBUG: built occ and neighcells for $(length(P)) pairs")
    println("DEBUG: picked $(picked)")
    println("Picked pairs and their letters:")
    for idx in picked
        println("Index: $idx, Letter: ", P[idx].t)
    end

    # If no pairs were generated, return early (no model to build)
    if length(P) == 0
        println("WARNING: no candidate pairs; skipping model build")
        return :NO_PAIRS, NaN, P, Int[]
    end

    if length(P) > 900 # TODO: comment out with HiGHS
        println("WARNING: too large number of pairs ($(length(P)));")
        return :NO_PAIRS, NaN, P, Int[]
    end

    # Sanity and preview: show greedy picks before building the model
    @assert all(p -> 1 ≤ p ≤ length(P), picked) "pick_greedy! returned out-of-range indices"
    if !isempty(picked)
        start_score = sum(P[i].weight for i in picked)
        println("Greedy preview (k=$(length(picked))) — starting sum = ", start_score)
        print_solution_grid(P, picked)
    end

    # Build and solve model inside try/catch to avoid throwing when solver/setup fails
    status = :ERROR
    obj = NaN
    chosen = Int[]
    try
        # model = Model(HiGHS.Optimizer)  # or GLPK.Optimizer / Cbc.Optimizer / Gurobi.Optimizer
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "time_limit", 120.0)  # 2-minute time limit
        # set_optimizer_attribute(model, "objective_bound", -174.0)
        # set_silent(model)

    @variable(model, y[1:length(P)], Bin)

    # Fix greedy-picked pairs (do not use start values)
    for p in picked
        @constraint(model, y[p] == 1)
    end

    # at most one pair per letter (disallow any skipping)
    for t in LETTERS
        ids = index[t]
        if !isempty(ids)
            @constraint(model, sum(y[i] for i in ids) == 1)
        end
    end

    # no overlap is already implied by the conflicts set built from overlap;
    # but adding cell-wise ≤1 can strengthen. Build once:
    # (Optional) Build cell-wise cap tightening:
    # ...

    # Enforce non-overlap by cell-wise constraints using occ (only occupied cells)
    # Deduplicate indices per cell, avoid shadowing JuMP variable `y` by using xi, yi.
    # for xi in 1:W, yi in 1:H
    #     ids = unique(occ[xi, yi])
    #     if length(ids) > 1
    #         @constraint(model, sum(y[i] for i in ids) ≤ 1)
    #     end
    # end

    # Grouped neighborhood constraints for every other cell (2,2), (2,4), ...
    # For each grid cell (i,j) with step 2, gather occ[] indices from the
    # 3x3 neighborhood centered at (i,j) and require at most one selected pair
    # among all pairs occupying those neighbor cells.
    
    for i in 2:2:W-1, j in 2:2:H-1
        ids = Int[]
        for dx in -1:1, dy in -1:1
            nx, ny = i + dx, j + dy
            if 1 ≤ nx ≤ W && 1 ≤ ny ≤ H
                append!(ids, occ[nx, ny])
            end
        end
        ids = unique(ids)
        if !isempty(ids)
            @constraint(model, sum(y[k] for k in ids) ≤ 1)
        end
    end

    # For each pair, build a compact constraint that forbids selecting that pair
    # together with any pair occupying its neighbor cells. This reduces the number
    # of constraints compared to enumerating all pairwise conflicts.
    for pid in eachindex(P)
        # gather unique neighbor pair indices from occ at the neighbor cells
        qs = Int[]
        for nb in neighcells[pid]
            nx, ny = nb
            append!(qs, occ[nx, ny])
        end
        qs = unique(qs)
        # remove self and same-letter pairs
        filter!(q -> q != pid && P[q].t != P[pid].t, qs)
        M = min(12, length(qs)) # 6 "should" be 12, but let's be real ain't nobody putting 12 Pentaminoes that close together
        if !isempty(qs)
            @constraint(model, sum(y[q] for q in qs) ≤ M * (1 - y[pid]))
        end
    end

        @objective(model, Max, sum(P[i].weight * y[i] for i in 1:length(P)))

        # Warm start using pairs that exactly match EXAMPLE_SOLUTION
        # if !isempty(WS)
        #     for i in 1:length(P)
        #         set_start_value(y[i], 0)
        #     end
        #     for pid in WS
        #         set_start_value(y[pid], 1)
        #     end
        #     println("Warm start with ", length(WS), " letters; nominal score = ", sum(P[i].weight for i in WS))
        # end

        if optimize
                
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
        end
    catch err
        # Log and return a safe failure result instead of throwing
        println("ERROR: building or solving model failed: ", err)
        return :ERROR, NaN, P, Int[]
    end
    return status, obj, P, chosen
end

# Simple runner when executing this file

function run(greedy_k::Int = 0, optimize::Bool = true, randomize::Bool = false, save_to_file::Bool = false)
    status, best, P, chosen = solve_pentomino(A, greedy_k, optimize, randomize)
    
    println("status: ", status)
    println("objective: ", best)
    println("selected pairs: ", length(chosen))
    
    if !isempty(chosen)
        if save_to_file
            file_name = "solution_$(best).txt"
            println("Saving solution to $file_name")
            open(file_name, "w") do io
                println(io, "status: ", status)
                println(io, "objective: ", best)
                println(io, "selected pairs: ", length(chosen))
                println(io, "solution grid (reduced):")
                print_solution_grid(P, chosen; io=io)
            end
        else
            println("solution grid (reduced):")
            print_solution_grid(P, chosen)
        end
    end
    # success only when the solver proved optimality
    success = (status == MOI.OPTIMAL)
    return success, best
end

function can_fill_holes(touch_counts::Matrix{Int}, cur_letter_idx::Int)
    # quick check: count empty cells not touching any placed letter
    bio::Matrix{Int} = zeros(Bool, W, H) # TODO: can reuse it with global cookie cntr

    can_place_cnt = 0

    for x in 1:2:W, y in 1:2:H
        if bio[x,y] == 0 && bio[x,y+1] == 0 && bio[x+1,y] == 0 && bio[x+1,y+1] == 0 &&
            touch_counts[x,y] == 0 && touch_counts[x,y+1] == 0 && touch_counts[x+1,y] == 0 && touch_counts[x+1,y+1] == 0
            
            start = (x,y)
            q = [start]
            q2 = [start]
            bio[x,y] = 1
            cnt = 1
            # bfs flood fill to find connected component size
            while !isempty(q)
                (cx,cy) = pop!(q)
                for (dx,dy) in [(-2,0),(2,0),(0,-2),(0,2)]
                    nx, ny = cx+dx, cy+dy
                    if 1 ≤ nx ≤ W && 1 ≤ ny ≤ H && bio[nx,ny] == 0 && touch_counts[nx,ny] == 0 && touch_counts[nx,ny+1] == 0 && touch_counts[nx+1,ny] == 0 && touch_counts[nx+1,ny+1] == 0
                        bio[nx,ny] = 1
                        push!(q, (nx,ny))
                        push!(q2, (nx,ny))
                        cnt += 1
                    end
                end
            end
            while !isempty(q2)
                (cx,cy) = pop!(q2)
                for dx in -1:2, dy in -1:2
                    nx, ny = cx+dx, cy+dy
                    if 1 ≤ nx ≤ W && 1 ≤ ny ≤ H && bio[nx,ny] == 0 && touch_counts[nx,ny] == 0
                        bio[nx,ny] = 1
                        cnt += 1
                    end
                end
            end
            pentamino_size = 48 # most shapes, 5*4 + 28 neighbors, we place smaller P (44) early on
            component_cnt = cnt ÷ pentamino_size
            can_place_cnt += component_cnt

            # if can_place_cnt >= length(LETTERS) - cur_letter_idx + 1
            #    return true
            # end
        end
    end

    # println("DEBUG: can_place_cnt = $can_place_cnt for cur_letter_idx = $cur_letter_idx")
    # for x in 1:W
    #     for y in 1:H
    #         if bio[x,y] == 0 && touch_counts[x,y] == 0
    #             print('.')
    #         elseif touch_counts[x,y] > 0
    #             print('X')
    #         elseif bio[x,y] > 0
    #             print('-')
    #         end
    #     end
    #     println()
    # end
    # println()

    return can_place_cnt >= length(LETTERS) - cur_letter_idx + 1
end

function estimate_score(pairs::Dict{Symbol, Vector{NamedTuple}})
    s = 0
    for t in LETTERS
        if isempty(pairs[t])
            return -Inf
        end
        s += pairs[t][1].weight
    end
    return s
end

const node_count = Ref(0)
const cur_letters::Array{Symbol} = LETTERS # Fixed letter order
const cutoff = 170

## Greedy DFS for recursive letter selection and model solve
function greedy_dfs(
    pairs::Dict{Symbol, Vector{NamedTuple}},
    cur_letter_idx::Int, 
    depth,
    cur_score::Int = 0,
    fixed_chosen::Vector{Symbol} = Symbol[],
    touch_counts::Matrix{Int} = zeros(Int, W, H),
    save_to_file::Bool = true,
)
    global node_count, cur_letters, cutoff

    best_possible_score = estimate_score(pairs)
    # println("DFS at depth $depth, letter index $cur_letter_idx, current score $cur_score, best possible score $best_possible_score")

    if best_possible_score <= cutoff
        return :NO_PAIRS, NaN, NamedTuple[], Int[]
    end
    
    branch_factor = 15 # how many top candidates to try per letter
    
    if depth <= 0 && !can_fill_holes(touch_counts, cur_letter_idx) # prune branches that dont have enough empty space
        return :NO_SOLUTION, NaN, [], Int[]
    end

    node_count[] += 1
    if node_count[] % 100 == 0
        println("DEBUG: DFS at depth $depth, node $(node_count[])")
    end

    if depth <= 0
        status = :ERROR
        obj = NaN
        chosen = Int[]
        # Build conflicts from the pruned pairs
        P, index, occ, neighcells, WS = build_conflicts(pairs)
        if length(P) == 0
            return :NO_PAIRS, NaN, P, Int[]
        end
        println("DFS at depth 0: building model with ", length(P), " pairs")
        try
            model = Model()
            @variable(model, y[1:length(P)], Bin)
            # Fix variables for letters already chosen in upper DFS levels
            # exactly one pair per every other letter
            for t in LETTERS
                ids = index[t]
                if t in fixed_chosen
                    push!(chosen, ids[1])
                    @constraint(model, y[ids[1]] == 1)
                else
                    @constraint(model, sum(y[i] for i in ids) == 1)
                end
            end
            # 3x3 neighborhood constraints
            for i in 2:2:W-1, j in 2:2:H-1
                ids = Int[]
                for dx in -1:1, dy in -1:1
                    nx, ny = i + dx, j + dy
                    if 1 ≤ nx ≤ W && 1 ≤ ny ≤ H
                        append!(ids, occ[nx, ny])
                    end
                end
                ids = unique(ids)
                if !isempty(ids)
                    @constraint(model, sum(y[k] for k in ids) ≤ 1)
                end
            end
            # Small-M constraints for pairs' neighbors
            for pid in eachindex(P)
                qs = Int[]
                for nb in neighcells[pid]
                    nx, ny = nb
                    append!(qs, occ[nx, ny])
                end
                qs = unique(qs)
                filter!(q -> q != pid && P[q].t != P[pid].t, qs)
                M = min(12, length(qs))
                if !isempty(qs)
                    @constraint(model, sum(y[q] for q in qs) ≤ M * (1 - y[pid]))
                end
            end
            @objective(model, Max, sum(P[i].weight * y[i] for i in 1:length(P)))

            nvar = JuMP.num_variables(model)
            n_le = MOI.get(model, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}}())
            n_eq = MOI.get(model, MOI.NumberOfConstraints{MOI.ScalarAffineFunction{Float64}, MOI.EqualTo{Float64}}())
            ncons = n_le + n_eq

            if max(nvar, ncons) > 1000
                set_optimizer(model, HiGHS.Optimizer)
            else
                set_optimizer(model, CPLEX.Optimizer)
            end

            # print_solution_grid(P, chosen)

            set_silent(model)
            optimize!(model)
            status = termination_status(model)
            try
                obj = objective_value(model)
                best = round(Int, obj)
                chosen = findall(i -> value(y[i]) > 0.5, 1:length(P))
                
                if save_to_file && length(chosen) == length(LETTERS)
                    println("DFS: found solution with objective ", obj, " and ", length(chosen), " pairs")
                    file_name = "solution_$(best).txt"
                    println("Saving solution to $file_name")
                    open(file_name, "w") do io
                        println(io, "status: ", status)
                        println(io, "objective: ", best)
                        println(io, "selected pairs: ", length(chosen))
                        println(io, "solution grid (reduced):")
                        print_solution_grid(P, chosen; io=io)
                    end
                end
            catch _err
                obj = NaN
                chosen = Int[]
                # println("DFS: failed to extract objective value")
            end
        catch err
            # println("ERROR: building or solving model failed") # err
            return :ERROR, NaN, P, Int[]
        end
        
        return status, obj, P, chosen
    end

    best_status = :NO_SOLUTION
    best_obj = -Inf
    best_P = []
    best_chosen = Int[]

    t = cur_letters[cur_letter_idx]
    candidates = pairs[t]
    
    if isempty(candidates)
        return :ERROR, NaN, [], Int[]
    end

    function pair_value(pr)
        return pr.weight / (length(pr.cellsU) + length(pr.cellsL))
    end

    # sorted = sort(candidates, by=pair_value, rev=true) # TODO: uncomment?
    sorted = candidates
    top_candidates = sorted[1:min(branch_factor, length(sorted))] # TODO: going to improve this?
    # println("  For letter ", t, ": trying ", length(top_candidates), " pairs")
    for (i, pr) in enumerate(top_candidates)
        # println("    Pair ", i, " for letter ", t, ": weight=", pr.weight, " cells=", length(union(pr.cellsU, pr.cellsL)))
        next_pairs = deepcopy(pairs)
        next_pairs[t] = [pr]
        next_score = cur_score + pr.weight
        # Track this fixed selection and update touch counts on inflated cells
        push!(fixed_chosen, t) 
        inflated_cells = inflate(pr)
        for (x,y) in inflated_cells
            # use [y,x] to match grid[y,x] convention
            if in_bounds((x,y))
                touch_counts[y, x] += 1
            end
        end
        blocked_region = inflate(pr)
        pruned_count = 0
        for other_t in LETTERS
            if other_t == t
                continue
            end
            if isempty(next_pairs[other_t])
                continue
            end
            before_count = length(next_pairs[other_t])
            keep = filter(other_pr -> !any(c -> c in blocked_region, union(other_pr.cellsU, other_pr.cellsL)), next_pairs[other_t])
            next_pairs[other_t] = keep
            pruned_count += before_count - length(keep)
        end
        # println("      Pruned ", pruned_count, " incompatible pairs")
        status, obj, P, chosen = greedy_dfs(next_pairs, cur_letter_idx + 1, depth - 1, next_score, fixed_chosen, touch_counts, save_to_file)
        # Backtrack: undo touch counts and fixed selection
        for (x,y) in inflated_cells
            if in_bounds((x,y))
                touch_counts[y, x] -= 1
            end
        end
        pop!(fixed_chosen)
        if isfinite(obj) && obj > best_obj
            best_status = status
            best_obj = obj
            best_P = P
            best_chosen = chosen
            # println("      New best at depth ", 5 - depth, ": ", obj)
        end
    end

    return best_status, best_obj, best_P, best_chosen
end

# Wrapper to run greedy_dfs and print results
function run_greedy_dfs(depth = 5, save_to_file = true)
    global node_count, cur_letters
    t0 = time_ns()
    node_count[] = 0
    pairs = deepcopy(PAIRS)
    cur_letters = [:W, :U, :I, :P, :L, :X, :N, :Z, :T, :F, :V, :Y]
    # Reusable state for DFS
    fixed_chosen = Vector{Symbol}()
    touch_counts = fill(0, H, W)
    status, obj, P, chosen = greedy_dfs(pairs, 1, depth, 0, fixed_chosen, touch_counts, save_to_file)
    println("========================================")
    println("Total DFS nodes visited: ", node_count[])
    println("Greedy DFS search results:")
    println("Status: ", status)
    println("Objective: ", obj)
    println("Selected pairs: ", length(chosen))
    if !isempty(chosen)
        println("Solution grid (reduced):")
        print_solution_grid(P, chosen)
    end
    dt_s = (time_ns() - t0) / 1e9
    println(@sprintf("run_greedy_dfs total time: %.2f s", dt_s))
end

# Keep running randomized greedy (k=5) until an optimal solution is proven or max_tries is reached
function run_many(greedy_k::Int = 5, max_tries::Int = 5000; optimize::Bool = true, randomize::Bool = true, stop_on_success::Bool = false)
    total_ns = 0.0
    successes = 0
    best = -Inf
    best_run = 0
    last_score = NaN

    for attempt in 1:max_tries
        print("Attempt #$attempt: ")
        t0 = time_ns()
        success, score = run(greedy_k, optimize, randomize, true)
        last_score = score
        dt_ns = (time_ns() - t0)
        total_ns += dt_ns
        avg_ms = (total_ns / attempt) / 1e6
        dt_ms = dt_ns / 1e6
        println(@sprintf("%.2f ms | aggregate avg %.2f ms", dt_ms, avg_ms))

        successes += success ? 1 : 0
        if isfinite(score) && (score > best)
            best = score
            best_run = attempt

            file_name = "best_many_" + string(best) + ".txt"
            println("New best score $best on attempt #$attempt; saving to $file_name")
        end

        if stop_on_success && success
            println("Succeeded on attempt #$attempt with score = ", score)
            println(@sprintf("Aggregate average after %d run(s): %.2f ms", attempt, avg_ms))
            return true, score, attempt
        end
    end

    final_avg_ms = (total_ns / max_tries) / 1e6
    if stop_on_success
        println("Did not reach proven optimality in ", max_tries, " attempts")
        println(@sprintf("Aggregate average after %d run(s): %.2f ms", max_tries, final_avg_ms))
        return false, last_score, max_tries
    else
        println(@sprintf("Summary: avg=%.2f ms over %d runs | successes=%d/%d | best=%s on run %d",
                         final_avg_ms, max_tries, successes, max_tries, string(best), best_run))
        return final_avg_ms, successes, best
    end
end

return 0