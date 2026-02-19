using Pkg
Pkg.activate(".")

using JuMP, Ipopt
using CSV, DataFrames, JSON, Revise

using Plots
using LaTeXStrings
using Measures
default(fontfamily="times", titlefontsize = 9, guidefontsize = 9)

# Set seed for reproducibility
using Random
seed = 1234
Random.seed!(seed)

# tickfontsize, legendfontsize

includet("../src/util.jl")
includet("../src/system_dyn_JuMP.jl")
includet("../src/models_JuMP.jl")
includet("../src/MPC_scenarios_JuMP.jl")
includet("../src/plots_JuMP.jl")
includet("../src/simulate_JuMP.jl")

# ----------- Load parameters --------

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

# Define parameter tuple
p = (NA_params = NA_params, ND_params = ND_params, s = s, V = V)

# Evaluate Nd for a range of Q values
Q_val = range(20,40, length=100)
Nd_vals = [evaluate_NN(ND_params, [2.63, 0.3, 0.3, Q])[1] for Q in Q_val]
plot(Q_val, Nd_vals, xlabel="Q", ylabel="Nd", title="Nd vs Q")
F_val = range(0.2, 0.45, length=100)
Nd_vals_F = [evaluate_NN(ND_params, [2.63, 0.3, F, 28])[1] for F in F_val]
plot(F_val, Nd_vals_F, xlabel="F", ylabel="Nd", title="Nd vs F", ylims = (300,800))

## Reference control
scen_funs = [define_scenario01A, define_scenario01B, define_scenario01C]

for scen_fun in scen_funs

    p_scen = scen_fun()
    model = generate_JuMP_model_ref(p_scen, p)
    optimize!(model)
    println(solution_summary(model))
    plt = plot_ref_MPC(model, p_scen, p)
    savefig(plt, "figures/MPC_" * p_scen.name * ".pdf")
    plt    
end

## Scenario parameters
# Scenario 02A - EMPC with optimized F and Q
p02A = define_scenario02A()

# Generate JuMP model for EMPC scenario 02A
model_2A = generate_JuMP_model_EMPC(p02A, p)

optimize!(model_2A)
println(solution_summary(model_2A))

plt = plot_EMPC(model_2A, p02A, p)
plt
# Save as PDF
savefig(plt, "figures/MPC_scenario02A.pdf")
plt

p02B = (; p02A..., lambda = [1e4, 1e2, 1e6]) # Keep control input and slack penalties, but remove price from objective
model_2B = generate_JuMP_model_EMPC_no_price(p02B, p)
optimize!(model_2B)
println(solution_summary(model_2B))
plt_2B = plot_EMPC(model_2B, p02B, p)
savefig(plt_2B, "figures/MPC_scenario02B.pdf")

p02C = (; p02A..., eta_min = 0.70) # Increase minimum capture efficiency requirement
model_2C = generate_JuMP_model_EMPC(p02C, p)
optimize!(model_2C)
println(solution_summary(model_2C))
plt_2C = plot_EMPC(model_2C, p02C, p)
savefig(plt_2C, "figures/MPC_scenario02C.pdf")  

p02D = (; p02C..., lambda = [1e4, 1e2, 1e6]) # Keep control input and slack penalties, but remove price from objective
model_2D = generate_JuMP_model_EMPC_no_price(p02D, p)
optimize!(model_2D)
println(solution_summary(model_2D))
plt_2D = plot_EMPC(model_2D, p02D, p)
savefig(plt_2D, "figures/MPC_scenario02D.pdf")