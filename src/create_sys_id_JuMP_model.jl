function create_JuMP_model_sys_id(datasets, NN_settings, x0)

    # Create JuMP model object
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 5)
    set_optimizer_attribute(model, "max_iter", 1_000)
    set_optimizer_attribute(model, "tol", 1e-8)

    # Storage for all NN parameter variables
    nn_params = Dict()

    for (nn_name, NN) in pairs(NN_settings)
        n_neurons    = NN.n_neurons
        n_hidden     = length(n_neurons)
        input_size   = length(NN.input.mu)
        output_size  = length(NN.output.mu)

        nn_params[nn_name] = Dict()
        n_neurons_prev = input_size

        layer_i = 0
        for i in 1:n_hidden
            layer_i += 1
            n_i = Int(n_neurons[i])

            W = @variable(model, [1:n_i, 1:n_neurons_prev], 
                        base_name = "W_$(nn_name)_$(layer_i)")
            b = @variable(model, [1:n_i],                   
                        base_name = "b_$(nn_name)_$(layer_i)")

            set_start_value.(W, 0.5 * randn(size(W)))
            set_start_value.(b, 0.5 * randn(size(b)))

            nn_params[nn_name][Symbol("W_$(layer_i)")] = W
            nn_params[nn_name][Symbol("b_$(layer_i)")] = b

            n_neurons_prev = n_i
        end
        layer_i += 1

        # Output layer
        W = @variable(model, [1:output_size, 1:n_neurons_prev], 
                        base_name = "W_$(nn_name)_$(layer_i)")
        b = @variable(model, [1:output_size],                   
                        base_name = "b_$(nn_name)_$(layer_i)")

        set_start_value.(W, 0.5 * randn(size(W)))
        set_start_value.(b, 0.5 * randn(size(b)))

        nn_params[nn_name][Symbol("W_$(layer_i)")] = W
        nn_params[nn_name][Symbol("b_$(layer_i)")] = b
    end

    # Parameter NamedTuple
    p = (NN = NN_settings, params = nn_params, s = (sa=0.0, sd=0.0, st=0.0, s_yga=0.0), V = (Va=0.17, Vd=0.17, Vt=0.17))



    # Model definition JuMP
    nd = length(datasets)
    ny = size(datasets[1].Y, 2)
    nx = length(x0)

    # --- Per-dataset state trajectories ---
    # X[d] is (nx x N_d) — each dataset has its own independent trajectory
    x_lb = [0.0, 0.0, 0.0] # Lower bounds for states
    x_ub = [5.0, 5.0, 5.0] # Upper bounds for states

    X = Vector{Matrix{VariableRef}}(undef, nd)
    for d in 1:nd
        N_d = datasets[d].N
        
        X[d] = @variable(model, [i=1:nx, k=1:N_d],
                    lower_bound = x_lb[i], upper_bound = x_ub[i],
                    base_name = "X_$(d)")
    end

    # Warm start states
    for d in 1:nd
        for i in 1:nx
            set_start_value.(X[d][i, :], x0[i])
        end
    end



    # --------------------------------------------------------
    # Euler collocation with sub-steps between data points
    # --------------------------------------------------------

    dt_euler = 10/3600

    for d in 1:nd
        N_d   = datasets[d].N
        t_d  = datasets[d].t          # time 
        u_d   = datasets[d].u           # (nu × N_d)
        dist_d = datasets[d].dist       # disturbances (nd_vars × N_d)
        dt_d = diff(t_d)                     # time intervals between data points




        # Loop over data intervals
        for k in 1:N_d-1

            # Number of Euler sub-steps per data interval
            n_steps_k = max(1, round(Int, dt_d[k] / dt_euler))
            dt_sub  = dt_d[k] / n_steps_k        # actual sub-step (adjusted to fit exactly)

            # Inputs/disturbances held constant over the interval (ZOH)
            u_k   = u_d[:, k]
            dist_k = dist_d[:, k]

            # Propagate from X[:, k] to X[:, k+1] using n_steps Euler steps
            x_cur = [X[d][i, k] for i in 1:nx]   # start of interval

            for step in 1:n_steps_k
                dx = system_equations_sys_id(x_cur, u_k, dist_k, p)

                if step < n_steps_k
                    # Intermediate states — introduce auxiliary variables
                    x_next = @variable(model, [1:nx],
                                base_name = "X_$(d)_$(k)_sub$(step)")
                    for i in 1:nx
                        @constraint(model, x_next[i] == x_cur[i] + dt_sub * dx[i])
                    end
                    x_cur = x_next
                else
                    # Final sub-step connects to the next data-point state variable
                    for i in 1:nx
                        @constraint(model, X[d][i, k+1] == x_cur[i] + dt_sub * dx[i])
                    end
                end
            end
        end
    end


    # --------------------------------------------------------
    # Objective: average MSE across datasets
    # --------------------------------------------------------
    total_mse = AffExpr(0.0)


    for d in 1:nd
        
        N_d  = datasets[d].N
        y_d  = datasets[d].Y
        ny   = size(y_d, 1)
        u_d  = datasets[d].u
        dist_d = datasets[d].dist

        dataset_sse = @expression(model,
            sum(
                begin
                    x_k   = X[d][:, k]
                    u_k   = u_d[:, k]
                    dist_k = dist_d[:, k]
                    y_hat = output_equations_sys_id(x_k, u_k, dist_k, p)
                    e = y_d[k] - y_hat
                    e^2
                end
                for k in 1:N_d
            )
        )

        total_mse += dataset_sse / (N_d * nd)
    end

    # Add regularization on NN parameters
    lambda_reg = 1e-3
    for (nn_name, layers) in nn_params
        for (layer_name, var) in layers
            total_mse += lambda_reg * sum(var[i]^2 for i in eachindex(var))
        end
    end

    # Add regularization to encourage steady state conditions at the first m time steps of each dataset
    m_steady = 5
    lambda_steady = 1e-2
    for d in 1:nd
        u_d = datasets[d].u
        dist_d = datasets[d].dist
        u_d_avg = mean(u_d[:, 1:m_steady], dims=2)
        dist_d_avg = mean(dist_d[:, 1:m_steady], dims=2)
        for k in 1:m_steady
            x_k   = X[d][:, k]
            u_k   = u_d_avg
            dist_k = dist_d_avg
            dx = system_equations_sys_id(x_k, u_k, dist_k, p)
            for i in 1:nx
                total_mse += lambda_steady * dx[i]^2
            end
        end
    end

    # Define objective
    @objective(model, Min, total_mse)

    return model, p, X
end
