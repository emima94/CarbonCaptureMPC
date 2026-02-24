# System identification

using JuMP
using Ipopt
using CSV, DataFrames, JSON, Revise, JLD2
using Plots
using LaTeXStrings
using Printf
using Random 
using Statistics

seed = 1234
Random.seed!(seed)

default(fontfamily="times", titlefontsize = 9, guidefontsize = 9)

includet("../src/util.jl")
includet("../src/system_dyn_JuMP.jl")
includet("../src/sys_id_util.jl")
includet("../src/create_sys_id_JuMP_model.jl")
includet("../src/simulate_JuMP.jl")

datasets = load_data_sys_id(;every_k = 6)

NN_settings = get_NN_settings()

x0 = [2.63, 0.3, 0.3] # Initial guess for states

# Training dataset idxs:
train_idxs = [2, 3, 4, 5, 6, 7, 8, 9, 10]
nd_train = length(train_idxs)

datasets_train = datasets[train_idxs]

model, p, X = create_JuMP_model_sys_id(datasets_train, NN_settings, x0);

#print(model)

# Optimize
optimize!(model)

# --------------------------------------------------------
# Extract results after optimization
# --------------------------------------------------------

println("Termination status : ", termination_status(model))
println("Objective (MSE)    : ", objective_value(model))

# --- Extract NN parameters as plain Float64 ---
nn_opt_params = Dict()

nn_jump_params = p.params

for (nn_name, layers) in nn_jump_params
    nn_opt_params[nn_name] = Dict()
    for (layer_name, var) in layers
        nn_opt_params[nn_name][layer_name] = value.(var)
    end
end

x0_list = Vector{Vector{Float64}}(undef, nd_train)
for d in 1:nd_train
    x0_list[d] = value.(X[d][:, 1])
end


# --- Pretty print ---
for (nn_name, layers) in nn_opt_params
    println("\n=== $nn_name ===")
    for (layer_name, val) in layers
        println("  $layer_name = ", val)
    end
end
for d in 1:nd_train
    println("\nInitial state for dataset $d: ", x0_list[d])
end

p_opt = (NN = NN_settings, params = nn_opt_params, s = p.s, V = p.V)

# Save results
@save "results/sys_id_results.jld2" p_opt x0_list






# # Evaluate on training data

# ffun = system_equations_sys_id

# function hfun02(x, u, dist, p)
#     return [output_equations_sys_id(x, u, dist, p)]
# end

# plt_list = []

# for d in 1:nd_train
#     x0_d = x0_list[d]
#     u_d = datasets_train[d].u
#     dist_d = datasets_train[d].dist
#     t_vec = datasets_train[d].t
#     dt_sim = 180/3600

#     X_sim, Z_sim = euler_simulation(ffun, hfun02, x0_d, u_d, dist_d, p_opt, t_vec, dt_sim)

#     plt = plot(t_vec[1:end-1], Z_sim[1,:], label="yga")
#     plot!(plt, t_vec, datasets_train[d].Y, label="yga_data", linestyle=:dash)
#     push!(plt_list, plt)
# end

# plot(plt_list..., layout=(2,5), size=(1200,600), title="System ID: YGA vs. time")

# # Evaluate Nd vs. Q and F
# Q_val = range(15,40, length=250)
# F_val = range(0.15, 0.45, length=250)
# ca = 2.60
# cd = 0.3
# Nd = zeros(length(F_val), length(Q_val))
# for (i, F) in enumerate(F_val)
#     for (j, Q) in enumerate(Q_val)
#         input_Nd = [ca, cd, F_val[i], Q_val[j]]
#         Nd[i,j] = nn_forward(input_Nd, NN_settings[:Nd], nn_opt_params[:Nd])[1] # Output is scalar
#     end
# end

# heatmap(F_val, Q_val, Nd', xlabel="F", ylabel="Q", title="Identified Nd vs. F and Q")

# # plot Nd vs. Q for fixed F
# i = [1, 50, 100, 150, 200, 250]
# plt_Q = plot()
# for idx in i
#     plot!(plt_Q, Q_val, Nd[idx, :], xlabel="Q", ylabel="Nd",label = "F = $(round(F_val[idx], digits=2))")
# end
# plt_Q

# plt_F = plot()
# for idx in i
#     plot!(plt_F, F_val, Nd[:, idx], xlabel="F", ylabel="Nd", label = "Q = $(round(Q_val[idx], digits=2))")
# end
# plt_F


