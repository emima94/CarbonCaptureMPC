# Load pkgs
using CSV, DataFrames, JSON

includet("../src/util.jl")
includet("../src/model.jl")

# Read parameter estimates from R
par_estimates = CSV.read("results/parameter_estimates_Tiller.csv", DataFrame)
par_fixed = CSV.read("results/parameter_fixed_Tiller.csv", DataFrame)
# Read JSON file
nn_settings = JSON.parsefile("results/nn_settings_Tiller.json")

# Convert to NamedTuple for easier access
par_estimates_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_estimates)]...)
par_fixed_nt = (; [Symbol(row.Column1) => row.x for row in eachrow(par_fixed)]...)

# Extract NN parameters
n_neurons_A = [1, 1]
n_neurons_D = [4, 8, 1]

NA_params = extract_NN_params(par_estimates, "a", n_neurons_A, nn_settings["Na"])
ND_params = extract_NN_params(par_estimates, "d", n_neurons_D, nn_settings["Nd"])

# Extract other parameters
s = (sa = par_estimates_nt.sa, sd = par_fixed_nt.sd, st = par_fixed_nt.st, s_yga = par_fixed_nt.s_yga)
V = (Va = par_fixed_nt.Va, Vd = par_fixed_nt.Vd, Vt = par_fixed_nt.Vt)

# Define input and disturbance tuples for testing
u = (F = 1.0, Q = 2.0)
d = (cgina = 1.0, Fgina = 2.0)
# Define parameter tuple
theta = (NA_params = NA_params, ND_params = ND_params, s = s, V = V)
# Complete parameter tuple
p = (u = u, d = d, theta = theta)

# Update p (by creating a new p, since p is immutable) (testing)
F_new = 7.0
u_new = (; p.u..., F = F_new)
p = (; p..., u = u_new)
p

# Test system equations
x = [2.5, 0.2, 0.2]
t = 0.0
dx = zeros(3)
system_equations!(dx, x, p, t)
dx
