using JuMP
using CPLEX  # or GLPK, HiGHS, Cbc, Gurobi

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

# helper for 8-neighborhood excluding the cell itself
function neighbors8((x,y))
    ((x-1,y-1),(x,y-1),(x+1,y-1),
        (x-1,y),(x+1,y),
        (x-1,y+1),(x,y+1),(x+1,y+1))
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
                        # Heuristic speed-up: skip negative placements
                        if total_sum >= 4 # 0 would be nice
                            push!(plist, (cellsU = cu, cellsL = cl, weight = total_sum, is_ws = matches_example))
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
    conflicts = Set{Tuple{Int,Int}}()

    # overlaps: any two different-letter pairs sharing a cell
    for x in 1:W, y in 1:H
        lst = occ[x,y]
        for i in 1:length(lst), j in (i+1):length(lst)
            p, q = lst[i], lst[j]
            if P[p].t != P[q].t
                push!(conflicts, (min(p,q), max(p,q)))
            end
        end
    end

    # neighbor touches: for each pair, look at neighbors of its occupied cells and conflict with pairs occupying those neighbor cells
    println("DEBUG: checking neighbor touches")

    # prepare per-pair neighbor-cell lists (this is what the caller expects)
    neighcells = [Vector{Tuple{Int,Int}}() for _ in 1:length(P)]

    for (pid, pr) in enumerate(P)
        # compute the outer-edge neighbor cells of the whole shape (unique)
        neighs = Set{Tuple{Int,Int}}()
        for c in pr.cells
            for nb in neighbors8(c)
                if in_bounds(nb)
                    push!(neighs, nb)
                end
            end
        end
        # remove cells that are actually occupied by the pair to get the "outer edge"
        for c in pr.cells
            if c in neighs
                delete!(neighs, c)
            end
        end

        # record neighbor cells for this pair (caller uses this)
        neighcells[pid] = collect(neighs)

        # iterate outer neighbors (unique) once and add conflicts with pairs occupying them
        for nb in neighs
            nx, ny = nb
            for q in occ[nx, ny]
                if pid < q && P[pid].t != P[q].t
                    push!(conflicts, (pid, q))
                end
            end
        end
    end

    return P, index, occ, neighcells, WS
end

# Build and solve the MILP
function solve_pentomino(A, run::Bool = false)
    ups   = enumerate_uppercase_placements()
    pairs = build_pairs(A, ups)
    println("DEBUG: built $(sum(length(pairs[t]) for t in LETTERS)) pairs")
    P, index, occ, neighcells, WS = build_conflicts(pairs)
    println("DEBUG: built occ and neighcells for $(length(P)) pairs")

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
        # --- Build pairwise conflicts from 3x3 pools and neighbor cells ---
        conflicts = Set{Tuple{Int,Int}}()
        # 3x3 pools centered on even-even cells
        for xi in 2:2:W-1, yi in 2:2:H-1
            pool = Int[]
            append!(pool, occ[xi, yi])
            for (nx,ny) in neighbors8((xi,yi))
                if 1 ≤ nx ≤ W && 1 ≤ ny ≤ H
                    append!(pool, occ[nx, ny])
                end
            end
            pool = unique(pool)
            if length(pool) > 1
                for i in 1:length(pool), j in (i+1):length(pool)
                    p = pool[i]; q = pool[j]
                    if P[p].t != P[q].t
                        push!(conflicts, (min(p,q), max(p,q)))
                    end
                end
            end
        end
        # Neighbor-touch pairs (small-M neighbors flattened to pairwise)
        for pid in 1:length(P)
            qs = Int[]
            for nb in neighcells[pid]
                nx, ny = nb
                append!(qs, occ[nx, ny])
            end
            qs = unique(qs)
            filter!(q -> q != pid && P[q].t != P[pid].t, qs)
            for q in qs
                push!(conflicts, (min(pid,q), max(pid,q)))
            end
        end

        # --- Write AMPL model (.mod) and data (.dat) using the new structure ---
        out_dir = normpath(joinpath(@__DIR__, ".."))
        mod_path = joinpath(out_dir, "pentomino.mod")
        dat_path = joinpath(out_dir, "pentomino.dat")

        open(mod_path, "w") do io
            println(io, "# Pentomino MILP (auto-generated; pairwise no-touch)\n")
            println(io, "set LETTERS;")
            println(io, "set Pairs;")
            println(io, "set PairsByLetter {l in LETTERS} within Pairs;\n")
            println(io, "set Conflicts within {i in Pairs, j in Pairs: i < j};")
            println(io, "set WarmStart within Pairs default {};\n")
            println(io, "param w {Pairs};\n")
            println(io, "var y {Pairs} binary;\n")
            println(io, "s.t. OnePerLetter {l in LETTERS}: sum {p in PairsByLetter[l]} y[p] = 1;")
            println(io, "s.t. NoTwo {(i,j) in Conflicts}: y[i] + y[j] <= 1;\n")
            println(io, "maximize OBJ: sum {p in Pairs} w[p] * y[p];\n")
        end

        open(dat_path, "w") do io
            println(io, "data;")
            println(io)

            # Letters
            print(io, "set LETTERS :=")
            for t in LETTERS
                print(io, " '", String(t), "'")
            end
            println(io, ";")
            println(io)

            # Pairs
            println(io, "set Pairs :=")
            for p in 1:length(P)
                print(io, p, (p % 20 == 0 ? "\n" : " "))
            end
            println(io, ";")
            println(io)

            # PairsByLetter
            for t in LETTERS
                print(io, "set PairsByLetter['", String(t), "'] :=")
                for pid in index[t]
                    print(io, " ", pid)
                end
                println(io, ";")
            end
            println(io)

            # Weights
            println(io, "param w :=")
            for p in 1:length(P)
                println(io, " ", p, " ", P[p].weight)
            end
            println(io, ";")
            println(io)

            # Pairwise conflicts
            println(io, "set Conflicts :=")
            ct = 0
            for (i,j) in sort(collect(conflicts))
                print(io, " (", i, ", ", j, ")")
                ct += 1
                if ct % 12 == 0
                    println(io)
                else
                    print(io, " ")
                end
            end
            println(io, ";")
            println(io)

            # No RowCuts/NeighPairs/M in pairwise format

            # Warm start indices
            print(io, "set WarmStart :=")
            for (k,p) in enumerate(WS)
                print(io, " ", p)
                if k % 30 == 0; print(io, "\n"); end
            end
            println(io, ";")
            println(io)
        end

    println("Exported AMPL files (pairwise):\n  MOD: ", mod_path, "\n  DAT: ", dat_path)
        # Continue to create a JuMP model and export LP/MPS as before
        model = Model(CPLEX.Optimizer)  # or GLPK.Optimizer / HiGHS.Optimizer / Cbc.Optimizer / Gurobi.Optimizer
        
        # CPLEX-specific attributes
        # set_optimizer_attribute(model, "CPX_PARAM_TILIM", 120.0)          # 2-minute time limit (seconds)
        set_optimizer_attribute(model, "CPXPARAM_Threads", max(1, Threads.nthreads()))
        # Optional: emphasize feasibility to accept/repair MIP starts faster
        # 0=balanced, 1=feasibility, 2=optimality, 3=best bound, 4=hidden feasibility
        set_optimizer_attribute(model, "CPX_PARAM_MIPEMPHASIS", 2)
        # Optional logging level (0=off, 1=milestones, 2=moderate)
        set_optimizer_attribute(model, "CPX_PARAM_MIPDISPLAY", 2)
        
        @variable(model, y[1:length(P)], Bin)
        # Assign stable, LP/MPS-safe variable names (y_1, y_2, ...)
        for i in 1:length(P)
            set_name(y[i], "y_" * string(i))
        end

    # at most one pair per letter (disallow any skipping)
    for t in LETTERS
        ids = index[t]
        if !isempty(ids)
            @constraint(model, sum(y[i] for i in ids) == 1)
        end
    end

    # Pairwise constraints
    for (i,j) in conflicts
        @constraint(model, y[i] + y[j] ≤ 1)
    end

    # no overlap is already implied by the conflicts set built from overlap;
    # but adding cell-wise ≤1 can strengthen. Build once:
    # (Optional) Build cell-wise cap tightening:
    # ...

    # No per-cell or small-M constraints; covered by pairwise Conflicts

    # For each pair, build a compact constraint that forbids selecting that pair
    # together with any pair occupying its neighbor cells. This reduces the number
    # of constraints compared to enumerating all pairwise conflicts.
    # (row-cuts removed)

        @objective(model, Max, sum(P[i].weight * y[i] for i in 1:length(P)))

        # Optionally record a warm start (note: LP/MPS won’t carry starts)
        if !isempty(WS)
            for i in 1:length(P)
                set_start_value(y[i], 0)
            end
            for pid in WS
                set_start_value(y[pid], 1)
            end
            println("Prepared warm start with ", length(WS), " letters; nominal score = ", sum(P[i].weight for i in WS))
        end

    # Export LP/MPS to same directory (already defined)
        #lp_path  = joinpath(out_dir, "pentomino_model.lp")
        #mps_path = joinpath(out_dir, "pentomino_model.mps")
        #JuMP.write_to_file(model, lp_path)
        #JuMP.write_to_file(model, mps_path)
        #println("Exported model to:\n  LP : " * lp_path * "\n  MPS: " * mps_path)

        if (run)
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

        status = :EXPORTED
        obj = NaN
        chosen = Int[]
    catch err
        # Log and return a safe failure result instead of throwing
        println("ERROR: building or solving model failed: ", err)
        return :ERROR, NaN, P, Int[]
    end
    return status, obj, P, chosen
end

# Simple runner when executing this file

function run(solve::Bool = false)
    status, best, _, _ = solve_pentomino(A, solve)
    if solve
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

    else
        println("status: ", status, " (model written to LP/MPS; not solved)")
    end
end

return 0