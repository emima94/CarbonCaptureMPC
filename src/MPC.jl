function solve_OCP(obj, U0, lb, ub, callback;
                   solver = OptimizationOptimJL.BFGS(), 
                   cons = nothing, 
                   lb_con = nothing, 
                   ub_con = nothing)
                   
    if cons !== nothing
        optf = OptimizationFunction(obj, DifferentiationInterface.SecondOrder(ADTypes.AutoForwardDiff(), ADTypes.AutoForwardDiff()), cons)
        opt_prob = OptimizationProblem(optf, U0, 0, lb=lb, ub=ub, lcons=lb_con, ucons=ub_con)
    else
        optf = OptimizationFunction(obj, ADTypes.AutoForwardDiff()) 
        opt_prob = OptimizationProblem(optf, U0, 0, lb=lb, ub=ub)
    end

    sol = Optimization.solve(opt_prob, solver, callback=callback, maxiters=500)
    return sol
end

