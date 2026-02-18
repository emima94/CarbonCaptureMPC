function simulate_control_intervals(U, D, prob, t_vec)
    T = eltype(U)

    x0 = prob.u0

    nx = length(x0)
    nz = 4

    X = zeros(T, nx, length(t_vec))
    Z = zeros(T, nz, length(t_vec)-1) 
    X[:,1] = x0
    
    for i in 1:(length(t_vec)-1)
        # Update p with current interval's inputs and disturbances
        u = (F = U[1,i], Q = U[2,i])
        d = (cgina = D[1,i], Fgina = D[2,i])
        p = (; prob.p..., u = u, d = d)

        # Remake problem with updated p
        prob_i = remake(prob; p = p, u0 = X[:,i], tspan = (t_vec[i], t_vec[i+1]))
        #solver = Tsit5()
        solver = RK4()
        ode_dt = (t_vec[i+1] - t_vec[i])/10
        sol_i = solve(prob_i, solver, dt = ode_dt, saveat=t_vec[i:i+1], adaptive = false)
        X[:,i+1] = sol_i.u[end]

        # Compute output
        Z[:,i] = output_equations(X[:,i], p, t_vec[i])

    end

    return X, Z
end

# Detailed version for plotting
function simulate_fine_grid(U, D, prob, t_vec, n_points_per_interval=10)
    T = eltype(U)
    nx, nz = 3, 4
    
    n_total = (length(t_vec) - 1) * n_points_per_interval + 1
    X = zeros(T, nx, n_total)
    Z = zeros(T, nz, n_total - 1)
    t_fine = zeros(T, n_total)
    
    X[:,1] = prob.u0
    t_fine[1] = t_vec[1]
    idx = 1
    
    for i in 1:(length(t_vec)-1)
        u = (F = U[1,i], Q = U[2,i])
        d = (cgina = D[1,i], Fgina = D[2,i])
        p = (; prob.p..., u = u, d = d)
        
        prob_i = remake(prob; p = p, u0 = X[:,idx], tspan = (t_vec[i], t_vec[i+1]))
        
        t_save = range(t_vec[i], t_vec[i+1], length=n_points_per_interval+1)
        ode_dt = (t_vec[i+1] - t_vec[i]) / n_points_per_interval
        sol_i = solve(prob_i, RK4(), dt=ode_dt, saveat=t_save, adaptive=false)
        
        for j in 2:length(sol_i.t)
            idx += 1
            X[:,idx] = sol_i.u[j]
            Z[:,idx-1] = output_equations(sol_i.u[j-1], p, sol_i.t[j-1])
            t_fine[idx] = sol_i.t[j]
        end
    end
    
    return X, Z, t_fine
end
