# Pentomino MILP (auto-generated; pairwise no-touch)

set LETTERS;
set Pairs;
set PairsByLetter {l in LETTERS} within Pairs;

set Conflicts within {i in Pairs, j in Pairs: i < j};
set WarmStart within Pairs default {};
set FixedOn within Pairs default {};

param w {Pairs};

var y {Pairs} binary;

s.t. OnePerLetter {l in LETTERS}: sum {p in PairsByLetter[l]} y[p] = 1;
s.t. NoTwo {(i,j) in Conflicts}: y[i] + y[j] <= 1;

s.t. FixOn {p in FixedOn}: y[p] = 1;

maximize OBJ: sum {p in Pairs} w[p] * y[p];

