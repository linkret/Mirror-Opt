
const N = 32

# Original 16x16 point matrix (kept in top-left); rest will be zero-padded to N x N
const _POINT_VALUE_MATRIX_ORIG = [
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

const EXAMPLE_SOLUTION_16 = [
    "iiiiiIIIII.v.T..";
    "...........v.TTT";
    "PP..ZZ...vvv.T..";
    "PP..Z....VVV.t..";
    "P..ZZ......V.ttt";
    "p..zz.llLL.V.t..";
    "pp..z..lL......n";
    "pp..zz.lL.U.U..n";
    ".......lL.UUU.nn";
    ".x..X.....uuu.n.";
    "xxxXXX.F..u.u.N.";
    ".x..X..FF.....NN";
    "......FF...yY..N";
    "..wW..ff..yyYY.N";
    ".wwWW..ff..yY...";
    "ww..WW.f...yY..."
]

# Build a 32x32 Char matrix by duplicating each character into a 2x2 block
const EXAMPLE_SOLUTION = begin
    rows16 = EXAMPLE_SOLUTION_16
    mat = Array{Char}(undef, 32, 32)
    for i in 1:16
        s = collect(rows16[i])
        # expand horizontally: each char twice
        outrow = Char[]
        for c in s
            push!(outrow, c)
            push!(outrow, c)
        end
        r1 = 2*i - 1
        r2 = 2*i
        for j in 1:32
            mat[r1, j] = outrow[j]
            mat[r2, j] = outrow[j]
        end
    end
    mat
end

const POINT_VALUE_MATRIX = begin
    m = zeros(Int, N, N)
    h0, w0 = size(_POINT_VALUE_MATRIX_ORIG)
    # place each original cell into the bottom-left cell of its 2x2 block
    for i in 1:h0, j in 1:w0
        newx = 2*(j-1) + 1   # left column of the 2x2 block
        newy = 2*(i-1) + 1   # bottom row of the 2x2 block
        #if newx <= N && newy <= N
        m[newy, newx] = _POINT_VALUE_MATRIX_ORIG[i, j]
        #end
    end
    m
end

const A = copy(POINT_VALUE_MATRIX)
const H, W = N, N

# Derived pentomino base shapes from the uppercase placements in EXAMPLE_SOLUTION.
# Shapes are returned as 5 relative (x,y) offsets with min x/y at 0.
const SHAPES_BASE = Dict(
    :F => [(0, 1), (1, 0), (1, 1), (1, 2), (2, 0)], 
    :V => [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)], 
    :X => [(0, 1), (1, 0), (1, 1), (1, 2), (2, 1)], 
    :Z => [(0, 0), (0, 1), (1, 1), (2, 1), (2, 2)], 
    :N => [(0, 1), (1, 1), (2, 0), (2, 1), (3, 0)], 
    :Y => [(0, 1), (1, 1), (2, 0), (2, 1), (3, 1)], 
    :T => [(0, 2), (1, 0), (1, 1), (1, 2), (2, 2)], 
    :P => [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)], 
    :I => [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)], 
    :U => [(0, 0), (0, 1), (1, 0), (2, 0), (2, 1)], 
    :W => [(0, 2), (1, 1), (1, 2), (2, 0), (2, 1)], 
    :L => [(0, 0), (0, 1), (1, 1), (2, 1), (3, 1)]
)

function _expand_shape(shape::Vector{Tuple{Int,Int}})
    pts = Tuple{Int,Int}[]
    for (x,y) in shape
        push!(pts, (2*x,   2*y))
        push!(pts, (2*x+1, 2*y))
        push!(pts, (2*x,   2*y+1))
        push!(pts, (2*x+1, 2*y+1))
    end
    # normalize so min x/y == 0
    xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
    dx, dy = minimum(xs), minimum(ys)
    sort!([(x-dx, y-dy) for (x,y) in pts])
end

const SHAPES = Dict(k => _expand_shape(v) for (k,v) in SHAPES_BASE)

const LETTERS = [:F, :I, :L, :N, :P, :T, :U, :V, :W, :X, :Y, :Z]

const DIRS = Dict(:U => (0,-1), :D => (0,1), :L => (-1,0), :R => (1,0))