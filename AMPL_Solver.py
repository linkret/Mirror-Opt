from amplpy import AMPL, ampl_notebook
import os

with open("ampl.lic", "r") as f:
    license_string = f.read().strip()

# Instantiate AMPL object and register magics with CPLEX module
ampl = ampl_notebook(
    modules=["cplex"],  # module to install
    license_uuid=license_string,
)

# ampl = AMPL() # Local AMPL instance

# ampl.eval("reset;") # Probably not needed

model_filename = "pentomino.mod"
data_filename = "pentomino.dat"
ampl.read(model_filename)
ampl.read(data_filename)

# Warm Start
# Set warm start: 0 everywhere, 1 on WarmStart
ampl.eval(r"""
let {p in Pairs} y[p] := 0;
let {p in WarmStart} y[p] := 1;
""")

# Solve with CPLEX
ampl.setOption("solver", "cplex")
ampl.setOption("cplex_options", "threads=4 mipdisplay=2") # mipemphasis=1 mipdisplay=2 nodefileind=2

ampl.solve()

# Inspect
ampl.eval("display solve_result, _solve_elapsed_time;")
ampl.eval("display {p in Pairs: y[p] > 0.5} y[p];")