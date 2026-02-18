using Pkg
Pkg.activate(".")

using Optimization, OptimizationNLopt, ADTypes, DifferentiationInterface, ForwardDiff
using OptimizationOptimJL
using OptimizationMOI, Ipopt 
using LinearAlgebra

# ---------- Problem setup ----------
N = 10               # horizon
x0 = 2.0
ymax = 3.0
Q = 1.0
R = 0.1
ρ = 100.0            # slack penalty

# Decision variables:
# z = [u[1:N]; s[1:N]]
n_u = N
n_s = N
n_z = n_u + n_s

# ---------- Dynamics & output ----------
function rollout(u)
    T = eltype(u)
    x = zeros(T, N+1)
    x[1] = x0
    for k = 1:N
        x[k+1] = x[k] + u[k]
    end
    y = x[2:end]
    return x, y
end

# ---------- Objective ----------
function cost(z, p)
    u = z[1:N]
    s = z[N+1:end]

    x, y = rollout(u)

    J = sum(Q*y.^2) + sum(R*u.^2) + ρ*sum(s.^2)
    return J
end

# ---------- Output-dependent constraint ----------
function constraint!(g, z, p)
    u = z[1:N]
    s = z[N+1:end]

    _, y = rollout(u)

    # y_k - ymax - s_k ≤ 0
    for k = 1:N
        g[k] = y[k] - ymax - s[k]
    end
end
# Total number of constraints
n_g = N

# Evaluate cost and constraints at initial guess
z0 = ones(n_z)*0.5
x, y = rollout(z0[1:N])
println("Initial x: ", x)
println("Initial y: ", y)
J0 = cost(z0, nothing)
g0 = zeros(n_g)
constraint!(g0, z0, nothing)
println("Initial cost: ", J0)
println("Initial constraints: ", g0)



# ---------- Callback for optimization progress ----------
# Callback
iter_count = Ref(0)
function callback(state, loss_val)
    iter_count[] += 1
    if iter_count[] % 10 == 0
        println("Iteration: $(iter_count[]), Loss: $(round(loss_val, digits=4))")
    end
    return iter_count[] >= 500
end

# ---------- Box constraints ----------
u_min = -1.0
u_max =  1.0
s_min =  0.0
s_max =  10.0

lb = [fill(u_min, N); fill(s_min, N)]
ub = [fill(u_max, N); fill(s_max, N)]

# ---------- Build Optimization problem ----------
z0 = ones(n_z)*0.5
# ---------- Build problem ----------
# Explicitly provide SecondOrder AD to silence the warning and get correct behaviour
ad = DifferentiationInterface.SecondOrder(
    ADTypes.AutoForwardDiff(),   # outer (for Hessian-vector products)
    ADTypes.AutoForwardDiff()    # inner
)

optf = OptimizationFunction(cost, ad; cons = constraint!)

optprob = OptimizationProblem(
    optf, z0, nothing;
    lb    = lb,
    ub    = ub,
    lcons = fill(-Inf, n_g),
    ucons = zeros(n_g)
)

# Ipopt is a proper interior-point NLP solver that uses second-order info
sol = solve(optprob, Ipopt.Optimizer(); callback = callback, maxiters = 500)


u_opt = sol.u[1:N]
s_opt = sol.u[N+1:end]

_, y_opt = rollout(u_opt)

println("Optimal y: ", y_opt)
println("Slack:     ", s_opt)
