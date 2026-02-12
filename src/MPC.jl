function solve_OCP(obj, U0, lb, ub)
    optf = OptimizationFunction(obj, ADTypes.AutoForwardDiff()) 
    opt_prob = Optimization.OptimizationProblem(optf, U0, 0, lb=lb, ub=ub)

    sol = Optimization.solve(opt_prob, OptimizationOptimJL.BFGS(), callback=callback, maxiters=500)
    return sol
end