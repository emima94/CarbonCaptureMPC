# System identification, minimizing RMSE between predicted and measured output (yga)
using Pkg
Pkg.activate(".")


using DifferentialEquations
import SciMLSensitivity as SMS
using Optimization, OptimizationOptimisers
using OptimizationOptimJL
using ForwardDiff
using DataInterpolations
using ComponentArrays
using Statistics
using Plots
using Printf
using CSV, DataFrames, Revise, JLD2
using LaTeXStrings
using Random 

seed = 1234
Random.seed!(seed)

default(fontfamily="times", titlefontsize = 9, guidefontsize = 9)

includet("../src/util_diffeq_sys_id.jl")

# Load data
datasets = load_data_sys_id(;every_k = 6)

# Training dataset idxs:
train_idxs = [2, 3, 4, 5, 6, 7, 8, 9, 10]
nd_train = length(train_idxs)
datasets_train = datasets[train_idxs]

# NN settings:
NN_settings = get_NN_settings()
p_fixed     = get_NN_normalization(NN_settings)

# NN parameters to be estimated:
Na_params = initialize_NN_params_from_settings(NN_settings.Na)
Nd_params = initialize_NN_params_from_settings(NN_settings.Nd)

p0 = (Na_params = Na_params, Nd_params = Nd_params)

# Initial guess for initial states
x0_guess = [[2.63; 0.3; 0.3] for _ in 1:length(datasets_train)]

# All parameters as ComponentArray:
θ0 = ComponentArray((p = p0, x0 = x0_guess))

# Create ODE problem for each dataset
prob_list = [make_ode_problem(ffun, ds, x0_guess[i], p0, p_fixed) for (i, ds) in enumerate(datasets_train)];

#sensealg = SMS.AutoZygote()
sensealg = SMS.InterpolatingAdjoint(autojacvec = SMS.ZygoteVJP())

### Optimize! ###
loss_history = Float64[]

# Optimizer settings
ad_type = Optimization.AutoForwardDiff()

# Optimization function
optf = Optimization.OptimizationFunction((θ, _) -> loss(θ, datasets_train, prob_list, sensealg, p_fixed), ad_type)

# Step 1: Adam optimizer
optimizer_adam = OptimizationOptimisers.ADAM(0.05)
prob_adam = Optimization.OptimizationProblem(optf, θ0)

println("\nStep 1: Optimizing with ADAM...")
result_adam = Optimization.solve(prob_adam, optimizer_adam;
    maxiters = 500,
    callback = callback_opt_02)

adam_iters = length(loss_history)

# Step 2: LBFGS optimizer
optimizer_lbfgs = OptimizationOptimJL.LBFGS(
    linesearch = LineSearches.BackTracking()
)
prob_lbfgs = Optimization.OptimizationProblem(optf, result_adam.u)

println("\nStep 2: LBFGS...")
result_lbfgs = Optimization.solve(prob_lbfgs, optimizer_lbfgs;
    maxiters = 100,
    callback = callback_opt_02)

sol = result_lbfgs.u
p_opt = NamedTuple(sol.p)
x0_opt = copy.(sol.x0)

# Save results
time_str = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
@save "results/sys_id_diffeq_MSE_results_$time_str.jld2" p_opt x0_opt


