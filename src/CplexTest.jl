using JuMP

if !isdefined(Main, :MIRROROPT_INCLUDED)
    include("MirrorOpt.jl")
    const MIRROROPT_INCLUDED = true
end

# Export pairwise (i+j <= 1) AMPL model using data built by MirrorOpt
function export_ampl(greedy_k::Int; outdir = normpath(joinpath(@__DIR__, "..")), modname = "pentomino.mod", datname = "pentomino.dat")
    # MirrorOpt provides the following interface:
    #   enumerate_uppercase_placements()
    #   build_pairs(A, ups) -> Dict{Symbol, Vector{NamedTuple}}
    #   build_conflicts(pairs) -> (P, index, occ, neighcells, WS)

    ups   = enumerate_uppercase_placements()
    pairs = build_pairs(A, ups)  # A is a global constant matrix from MirrorOpt
    # Greedy prune incompatible pairs and select k fixed ones
    picked = pick_greedy!(greedy_k, pairs)
    # Now pairs is mutated; build flattened structures
    P, index, occ, neighcells, WS = build_conflicts(pairs)

    println("DEBUG: picked $(picked)")
    println("Picked pairs and their letters:")
    for idx in picked
        println("Index: $idx, Letter: ", P[idx].t)
    end

    # Sanity and preview: show greedy picks before building the model
    @assert all(p -> 1 ≤ p ≤ length(P), picked) "pick_greedy! returned out-of-range indices"
    if !isempty(picked)
        start_score = sum(P[i].weight for i in picked)
        println("Greedy preview (k=$(length(picked))) — starting sum = ", start_score)
        print_solution_grid(P, picked)
    end

    # If no pairs were generated, return early (no model to build)
    if length(P) == 0
        println("WARNING: no candidate pairs; skipping model build")
        return :NO_PAIRS, NaN, P, Int[]
    end

    # Build pairwise conflicts from 3x3 pools and neighbor cells (as in MirrorOpt's recent design)
    conflicts = Set{Tuple{Int,Int}}()
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

    mod_path = joinpath(outdir, modname)
    dat_path = joinpath(outdir, datname)

    open(mod_path, "w") do io
        println(io, "# Pentomino MILP (auto-generated; pairwise no-touch)\n")
        println(io, "set LETTERS;")
        println(io, "set Pairs;")
        println(io, "set PairsByLetter {l in LETTERS} within Pairs;\n")
        println(io, "set Conflicts within {i in Pairs, j in Pairs: i < j};")
        println(io, "set WarmStart within Pairs default {};")
        println(io, "set FixedOn within Pairs default {};\n")
        println(io, "param w {Pairs};\n")
        println(io, "var y {Pairs} binary;\n")
        println(io, "s.t. OnePerLetter {l in LETTERS}: sum {p in PairsByLetter[l]} y[p] = 1;")
        println(io, "s.t. NoTwo {(i,j) in Conflicts}: y[i] + y[j] <= 1;\n")
        println(io, "s.t. FixOn {p in FixedOn}: y[p] = 1;\n")
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

        # Warm start indices
        print(io, "set WarmStart :=")
        for (k,p) in enumerate(WS)
            print(io, " ", p)
            if k % 30 == 0; print(io, "\n"); end
        end
        println(io, ";")
        println(io)

        # Greedy fixed-on indices
        print(io, "set FixedOn :=")
        for (kk,p) in enumerate(picked)
            print(io, " ", p)
            if kk % 30 == 0; print(io, "\n"); end
        end
        println(io, ";")
        println(io)
    end

    println("Exported AMPL files (pairwise):\n  MOD: ", mod_path, "\n  DAT: ", dat_path)
    return mod_path, dat_path
end

return 0