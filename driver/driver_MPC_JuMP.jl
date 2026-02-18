using Pkg
Pkg.activate(".")

using JuMP, Ipopt
using CSV, DataFrames, JSON, Revise

using Plots
using LaTeXStrings
default(fontfamily="times", titlefontsize = 12)

includet("../src/util.jl")
includet("../src/system_dyn_JuMP.jl")

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



# ---------- Problem setup ----------
t0 = 0
tf = 6
dt = 0.25
N = Int((tf - t0) / dt)  # number of intervals
t = range(t0, tf, length=N+1)

Nx_per_interval = 10
Nx = N * Nx_per_interval + 1   # total number of state nodes
dt_sim = dt / Nx_per_interval
t_sim = range(t0, tf, length=Nx)

# What control interval does each state belong to?
interval_of = repeat(1:N, inner=Nx_per_interval)
push!(interval_of, N)  # last state belongs to last interval

# Define scenario
cgina = repeat([4.4], inner=N) # Disturbance sequence (cgina)
Fgina = repeat([180], inner=N) # Disturbance sequence (Fgina)
D = hcat(cgina, Fgina)'

# Previous control input (will be updated in closed-loop)
F_prev = 0.3
Q_prev = 28.0

# Minimum capture efficiency
eta_min = 0.90

# Electricity and CO2 prices
p_el = repeat([0.2, 1.0], inner=Int(N/2)) # EUR/kWh
p_CO2 = 0.07  #EUR/kg

# MPC tuning Parameters
lambda_F = 1e3
lambda_Q = 1e1
lambda_s = 1e6

# ODE rhs:
ffun = system_equations_w_cum_CO2

# ----- Model -------
model = Model(Ipopt.Optimizer)
set_attribute(model, "print_level", 0)

# Initial state
x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state
nx = length(x0)

# States (fine grid, Nx nodes)
@variable(model, x[1:nx, 0:(Nx-1)]) # 5 states: ca, cd, ct, Ma, Mina
# Set initial guess for states (constant at initial state)
X_guess = repeat(x0, 1, Nx) # Initial guess for states (constant at initial state)
set_start_value.(x, X_guess)

# Controls (coarse grid, N intervals)
# bounds
F_min, F_max = 0.22, 0.4
Q_min, Q_max = 22.0, 35.0
s_min = 0.0
@variable(model, F[1:N], lower_bound=F_min, upper_bound=F_max)
@variable(model, Q[1:N], lower_bound=Q_min, upper_bound=Q_max)
@variable(model, s, lower_bound=s_min) # Slack variable for constraint violation

# Initial guess for controls
F_guess = fill(0.3, N)
Q_guess = fill(28.0, N)
set_start_value.(F, F_guess)
set_start_value.(Q, Q_guess)
set_start_value(s, 0.0)

# constraints
@constraint(model, initial_state[i = 1:nx],  x[i,0] == x0[i]) # Initial condition

# dynamic constraints using Euler collocation
for i in 1:Nx-1
    # Determine which control interval this state belongs to
    k = interval_of[i] # control interval index for state node i-1
    # Evaluate dynamics at current state and control
    u_k = [F[k], Q[k]]
    d_k = D[:,k]
    f = ffun(x[:,i-1], u_k, d_k, p)
    for s in 1:nx 
        @constraint(model, x[s,i] == x[s,i-1] + dt_sim * f[s])
    end
end

# Output constraint at final time step (Ensure capture efficiency meets minimum requirement at average over the horizon)
@constraint(model, x[4,end] - x[5,end] * eta_min + s >= 0); # Ma(tf) - Mina(tf)*eta_min/100 + s >= 0



# Objective: 
# Assuming a constant CO2 price 
@objective(model, Min,
    sum(p_el[k] * Q[k] * dt for k = 1:N) 
        + p_CO2 * x[4,end] * dt * N 
        + lambda_F * ((F[1] - F_prev)^2 + sum((F[k] - F[k-1])^2 for k = 2:N)) 
        + lambda_Q * ((Q[1] - Q_prev)^2 + sum((Q[k] - Q[k-1])^2 for k = 2:N)) 
        + lambda_s * s^2
)


optimize!(model)

println(solution_summary(model))

termination_status(model)
primal_status(model)
objective_value(model)
# Extract values

x_opt = value.(x)
F_opt = value.(F)
Q_opt = value.(Q)
s_opt = value(s)

hfun = output_equations


function euler_simulation(ffun, hfun, x0, U, D, p, t_vec)
    nx = length(x0)
    X = zeros(nx, length(t_vec))
    X[:,1] = x0
    
    z0 = hfun(x0, U[:,1], D[:,1], p)
    nz = length(z0)
    Z = zeros(nz, length(t_vec)-1)
    
    for i in 1:(length(t_vec)-1)
        f = ffun(X[:,i], U[:,i], D[:,i], p)
        X[:,i+1] = X[:,i] + dt_sim * f
        Z[:,i] = hfun(X[:,i], U[:,i], D[:,i], p)
    end
    return X, Z
end

U_sim = repeat(hcat(F_opt, Q_opt)', inner=(1,Nx_per_interval))
Dsim = repeat(D, inner=(1,Nx_per_interval))

X, Z = euler_simulation(ffun, hfun, x0, U_sim, Dsim, p, t_sim)

# Compute cost and revenues

cost = p_el .* Q_opt .* dt
revenue = sum(reshape(p_CO2 * diff(X[4,:]) * dt_sim, Nx_per_interval, N), dims=1)
profit = revenue' - cost

# Plotting
plt_cap_eff = plot(t_sim[2:end], Z[3,:], seriestype=:steppre, linestyle=:dot, label=L"y^g_{a,ref}", title="Output", xlabel="", ylabel="Cap. eff. [%]", ylims=(20,110))

plt2 = plot(t_sim[2:end], Z[1:2,:]', labels=["yga" "ygina"], title="Output Trajectories with Optimized F", xlabel="Time (hr)", ylabel="Concentration (vol%)", ylims = (-2,14))
plt3 = plot(t_sim[2:end], Z[4,:], labels=["Na"],  title = "Capture rate", xlabel="Time (hr)", ylabel="Capture Rate (kg/hr)", ylims=(0,50))

plt_x = plot(t_sim, X[1:3,:]', labels=[L"c_a" L"c_d" L"c_t"], title="States", xlabel="", ylabel="CO2 liq. conc. [kmol/m3]", ylims = (-1,5))
plt_F = plot(t, vcat(F_opt[1], F_opt), seriestype = :steppre, title="Optimized Control Input", xlabel="Time (hr)", ylabel="Flow Rate", ylims=(0.2,0.45), label=L"F")
plot!(plt_F, t, fill(F_min, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
plot!(plt_F, t, fill(F_max, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")

plt_Q = plot(t, vcat(Q_opt[1], Q_opt), seriestype = :steppre, label=L"Q", title="", xlabel="", ylabel="Reboiler duty [kW]", ylims = (20,37))
plot!(plt_Q, t, fill(Q_min, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
plot!(plt_Q, t, fill(Q_max, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")

plt_cgina = plot(t_sim, vcat(Dsim[1,1], Dsim[1,:]), seriestype = :steppre, label=L"c_{in,a}^g", title="Disturbances", xlabel="", ylabel="Flue gas conc. [mol/m3]", ylims=(4,5))
plt_Fgina = plot(t_sim, vcat(Dsim[2,1], Dsim[2,:]), seriestype = :steppre, label=L"F_{in,a}^g", title="", xlabel="", ylabel="Liq. flow rate [m3/h]", ylims=(140,220))

plt_p_el = plot(t, vcat(p_el[1], p_el) * 1000, seriestype=:steppre, label=L"p_{el}", title="Electricity Price", xlabel="", ylabel="Price [EUR/MWh]", ylims=(0,1200))
plt_profit = plot(t, vcat(profit[1], profit), seriestype=:steppre, label=L"Profit", title="Profit over Time", xlabel="Time (hr)", ylabel="Profit [EUR]", ylims=(-10,10))

plt_cap_cum = plot(t_sim, X[4,:], seriestype=:steppre, label=L"M_a", title="Cumulative CO2 Captured", xlabel="Time (hr)", ylabel="Mass of CO2 Captured [kg]", ylims=(0,300))
hline!([eta_min * X[5,end]], linestyle=:dash, color=:red, label="Min. CO2 to capture")

scale = 900.0
H2W_ratio = 0.8

plt_out = plot(plt_cap_eff, plt_F, plt_cgina, plt_x, plt_Q, plt_Fgina, plt_p_el, plt_profit, plt_cap_cum,
    layout = (3,3), link=:x, 
    size = (scale, scale * H2W_ratio)
)
