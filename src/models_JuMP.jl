function generate_JuMP_model_ref(p_scenario, p)

    # Extract parameters and functions
    cap_eff_ref = p_scenario.cap_eff_ref
    lambda_F = p_scenario.lambda
    x0 = p_scenario.x0
    F_min = p_scenario.lb
    F_max = p_scenario.ub
    N = p_scenario.N
    Nx = p_scenario.Nx
    interval_of = p_scenario.interval_of
    D = p_scenario.D
    dt = p_scenario.dt
    dt_sim = p_scenario.dt_sim
    ffun = p_scenario.ffun
    F_prev = p_scenario.F_prev
    Q = p_scenario.Q
    hfun = p_scenario.hfun
    ode_solver = p_scenario.ODE_solver

    # ----- Model -------
    model = Model(Ipopt.Optimizer)
    set_attribute(model, "print_level", 5)

    # Initial state
    nx = length(x0)
    nz = length(hfun(x0, fill(0.0, 2), D[:,1], p))

    # States (fine grid, Nx nodes)
    @variable(model, x[1:nx, 0:(Nx-1)]) # 5 states: ca, cd, ct, Ma, Mina
    # Set initial guess for states (constant at initial state)
    X_guess = repeat(x0, 1, Nx) # Initial guess for states (constant at initial state)
    set_start_value.(x, X_guess)



    # Controls (coarse grid, N intervals)
    @variable(model, F[1:N], lower_bound=F_min, upper_bound=F_max)

    # Initial guess for controls
    F_guess = fill(0.3, N)
    set_start_value.(F, F_guess)

    # constraints
    @constraint(model, initial_state[i = 1:nx],  x[i,0] == x0[i]) # Initial condition

    # ODE solver
    if ode_solver == "euler"
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
    elseif ode_solver == "rk4"
        # Dynamic constraints using RK4 collocation with intermediate variables
           # # Intermediate RK4 variables for collocation
        @variable(model, _k1[1:nx, 1:(Nx-1)])
        @variable(model, _k2[1:nx, 1:(Nx-1)])
        @variable(model, _k3[1:nx, 1:(Nx-1)])
        @variable(model, _k4[1:nx, 1:(Nx-1)])
        set_start_value.(_k1, X_guess[:,1:end-1])
        set_start_value.(_k2, X_guess[:,1:end-1])
        set_start_value.(_k3, X_guess[:,1:end-1])
        set_start_value.(_k4, X_guess[:,1:end-1])

        # Dynamic constraints using RK4 collocation with intermediate variables
        for i in 1:Nx-1
            k = interval_of[i] # control interval index for state node i-1
            u_k = [F[k], Q[k]]
            d_k = D[:,k]
            x_i = x[:,i-1]

            @constraint(model, _k1[:,i] == ffun(x_i, u_k, d_k, p))
            @constraint(model, _k2[:,i] == ffun(x_i + dt_sim/2.0 * _k1[:,i], u_k, d_k, p))
            @constraint(model, _k3[:,i] == ffun(x_i + dt_sim/2.0 * _k2[:,i], u_k, d_k, p))
            @constraint(model, _k4[:,i] == ffun(x_i + dt_sim * _k3[:,i], u_k, d_k, p))
            @constraint(model, x[:,i] == x_i + dt_sim/6.0 * (_k1[:,i] + 2*_k2[:,i] + 2*_k3[:,i] + _k4[:,i]))
        end
    else
        error("Unsupported ODE solver: $ode_solver")
    end


    # Get capture efficiency at each time step
    @expression(model, cap_eff[j=1:Nx-1], 
                    hfun(x[:,j], 
                        [F[interval_of[j]], Q[interval_of[j]]], 
                        D[:,interval_of[j]], 
                        p)[3]
    )

    # Objective: 
    # Reference control of capture efficiency.
    @objective(model, Min,
        sum((cap_eff[j] - cap_eff_ref[j])^2 for j in 1:Nx-1) * dt_sim
            + lambda_F * ((F[1] - F_prev)^2 + sum((F[k] - F[k-1])^2 for k = 2:N)) 
    )

    return model

end

function generate_JuMP_model_EMPC(p_scenario, p)


    # Extract parameters and functions
    p_el = p_scenario.p_el
    p_CO2 = p_scenario.p_CO2
    lambda = p_scenario.lambda
    eta_min = p_scenario.eta_min
    x0 = p_scenario.x0
    lb = p_scenario.lb
    ub = p_scenario.ub
    N = p_scenario.N
    Nx = p_scenario.Nx
    interval_of = p_scenario.interval_of
    D = p_scenario.D
    dt = p_scenario.dt
    dt_sim = p_scenario.dt_sim
    t = p_scenario.t
    t_sim = p_scenario.t_sim
    ffun = p_scenario.ffun
    hfun = p_scenario.hfun
    F_prev = p_scenario.F_prev
    Q_prev = p_scenario.Q_prev
    
    lambda_F, lambda_Q, lambda_s = lambda

    # ----- Model -------
    model = Model(Ipopt.Optimizer)
    set_attribute(model, "print_level", 5)

    # Initial state
    
    nx = length(x0)

    # States (fine grid, Nx nodes)
    @variable(model, x[1:nx, 0:(Nx-1)]) # 5 states: ca, cd, ct, Ma, Mina
    # Set initial guess for states (constant at initial state)
    X_guess = repeat(x0, 1, Nx) # Initial guess for states (constant at initial state)
    set_start_value.(x, X_guess)

    # Controls (coarse grid, N intervals)
    # bounds
    F_min, Q_min, s_min = lb
    F_max, Q_max, s_max = ub

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
            + p_CO2 * x[4,end] 
            + lambda_F * ((F[1] - F_prev)^2 + sum((F[k] - F[k-1])^2 for k = 2:N)) 
            + lambda_Q * ((Q[1] - Q_prev)^2 + sum((Q[k] - Q[k-1])^2 for k = 2:N)) 
            + lambda_s * s^2
    )

    return model

end

function generate_JuMP_model_EMPC_no_price(p_scenario, p)
    

    # Extract parameters and functions
    p_el = p_scenario.p_el
    p_CO2 = p_scenario.p_CO2
    lambda = p_scenario.lambda
    eta_min = p_scenario.eta_min
    x0 = p_scenario.x0
    lb = p_scenario.lb
    ub = p_scenario.ub
    N = p_scenario.N
    Nx = p_scenario.Nx
    interval_of = p_scenario.interval_of
    D = p_scenario.D
    dt = p_scenario.dt
    dt_sim = p_scenario.dt_sim
    t = p_scenario.t
    t_sim = p_scenario.t_sim
    ffun = p_scenario.ffun
    hfun = p_scenario.hfun
    F_prev = p_scenario.F_prev
    Q_prev = p_scenario.Q_prev
    
    lambda_F, lambda_Q, lambda_s = lambda

    model = generate_JuMP_model_EMPC(p_scenario, p)

    Q = model[:Q]
    F = model[:F]
    s = model[:s]

    # Rewrite objective to only minimize electricity consumption without considering price
    @objective(model, Min,
        sum(Q[k]^2 for k = 1:N) * dt
            + lambda_F * ((F[1] - F_prev)^2 + sum((F[k] - F[k-1])^2 for k = 2:N)) 
            + lambda_Q * ((Q[1] - Q_prev)^2 + sum((Q[k] - Q[k-1])^2 for k = 2:N)) 
            + lambda_s * s^2
    )

    return model
end


