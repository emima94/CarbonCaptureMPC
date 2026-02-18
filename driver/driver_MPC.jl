# Load pkgs
using Pkg

Pkg.activate(".")

using CSV, DataFrames, JSON, Revise
using DifferentialEquations, DataInterpolations
using Ipopt
using DifferentiationInterface, ADTypes
using Optimization, OptimizationOptimJL
using Plots
using LaTeXStrings
default(fontfamily="times", titlefontsize = 12)

using ForwardDiff, ADTypes

includet("../src/util.jl")
includet("../src/model.jl")
includet("../src/simulate.jl")
includet("../src/loss_functions.jl")
includet("../src/MPC.jl")
includet("../src/MPC_scenarios.jl")
includet("../src/plots.jl")

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
out_allocate = (u_out = zeros(2), d_out = zeros(2)) # Placeholder for output allocation parameters, if needed

# Define parameter tuple
theta = (NA_params = NA_params, ND_params = ND_params, s = s, V = V, out_allocate = out_allocate)

# Define p tuple for ODE problem (will be updated in simulate function)
p = (u = (F = NaN, Q = NaN), d = (cgina = NaN, Fgina = NaN), theta = theta)


## Cheack out ND vs Q_val
Q_val = range(20,40, length=100)
Nd_vals = [evaluate_NN(ND_params, [2.63, 0.3, 0.3, Q])[1] for Q in Q_val]
plot(Q_val, Nd_vals, xlabel="Q", ylabel="Nd", title="Nd vs Q")
F_val = range(0.2, 0.45, length=100)
Nd_vals_F = [evaluate_NN(ND_params, [2.63, 0.3, F, 28])[1] for F in F_val]
plot(F_val, Nd_vals_F, xlabel="F", ylabel="Nd", title="Nd vs F", ylims = (300,800))

##

# Simulation Time
t0 = 0.0
tf = 6.0
dt = 15/60   # Lenght of each control interval (hr)
N = Int((tf - t0) / dt) # Number of control intervals
t_vec = range(t0, tf, length=N+1)


## Ref. control scenarios ##

# Scenario functions, ref. control
scen_funcs = [define_scenario01A, define_scenario01B, define_scenario01C]
scen_fig_names = ["scenario01A", "scenario01B", "scenario01C"]



# for i in eachindex(scen_funcs)
# #i = 3
#     def_scen_func = scen_funcs[i]
#     scen_fig_name = scen_fig_names[i]
#     ### Scenario definition ###
#     obj_scenario, F0, lb, ub, cons, lb_con, ub_con, scen_data = def_scen_func(t0, tf, p, N);

#     # Solve OCP
#     solver = OptimizationOptimJL.BFGS()
#     sol = solve_OCP(obj_scenario, F0, lb, ub, callback, solver = solver, cons = cons, lb_con = lb_con, ub_con = ub_con)

#     # Simulate with optimized F values
#     Usim = hcat(sol.u, scen_data.U0[2,:])' # Combine optimized F with original Q
#     Dsim = scen_data.D0

#     Xsim, Zsim, t_fine_sim = simulate_fine_grid(Usim, Dsim, scen_data.prob, scen_data.t_vec, 50)

#     plt = plot_results(Usim, Dsim, Xsim, Zsim, t_fine_sim, scen_data)

#     savefig(plt, "figures/MPC_" * scen_fig_name * ".pdf")
# end

## Economic MPC scenarios ###

obj_scenario, U0, lb, ub, cons, lb_con, ub_con, scenario_data = define_scenario02A(t0, tf, p, N) 

optf = OptimizationFunction(obj_scenario, DifferentiationInterface.SecondOrder(ADTypes.AutoForwardDiff(), ADTypes.AutoForwardDiff()), cons = cons)

opt_prob = OptimizationProblem(optf, U0, 0, lb=lb, ub=ub, lcons=lb_con, ucons=ub_con)

sol = solve(opt_prob, Ipopt.Optimizer())