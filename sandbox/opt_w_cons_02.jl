using Pkg
Pkg.activate(".")

using OptimizationBase, OptimizationMOI, OptimizationOptimJL, Ipopt
using ForwardDiff
using DifferentiationInterface, ADTypes

rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2

x0 = zeros(2)
_p = [1.0, 1.0]

cons(res, x, p) = (res .= [x[1]^2 + x[2]^2, x[1] * x[2]])

optprob = OptimizationFunction(rosenbrock, DifferentiationInterface.SecondOrder(ADTypes.AutoForwardDiff(), ADTypes.AutoForwardDiff()), cons = cons)
prob = OptimizationProblem(optprob, x0, _p, lcons = [-Inf, -1.0], ucons = [0.8, 2.0])
sol = solve(prob, IPNewton())

res = zeros(2)
cons(res, sol.u, _p)
res

sol = solve(prob, Ipopt.Optimizer())

