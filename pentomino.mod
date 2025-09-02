# Pentomino MILP (auto-generated)

set LETTERS;
set Pairs;
set PairsByLetter {l in LETTERS} within Pairs;
set Conflicts within Pairs cross Pairs;
set WarmStart within Pairs;  # optional, for MIP starts

param w {Pairs};

var y {Pairs} binary;

s.t. OnePerLetter {l in LETTERS}: sum {p in PairsByLetter[l]} y[p] = 1;
s.t. NoConflict {(i,j) in Conflicts}: y[i] + y[j] <= 1;

maximize OBJ: sum {p in Pairs} w[p] * y[p];

# Use WarmStart in a .run script, e.g.:
# let {p in Pairs} y[p] := if p in WarmStart then 1 else 0;

