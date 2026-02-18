
function define_scenario01A(t0, tf, p, N) 
    ### Scenario 1 ###

    t_vec = range(t0, tf, length=N+1)

    # Define input and disturbance sequence for scenario
    cgina_val = [4.4] 
    Fgina_val = [160, 180]
    F_val = [0.26]
    Q_val = [28]
    cap_eff_ref_val = [90]
    cap_eff_ref = repeat(cap_eff_ref_val, inner=Int(N/length(cap_eff_ref_val)))
    # MPC tuning parameters
    lambda_dU = 1e5

    # Initial state
    x0 = [2.63; 0.3; 0.3]

    U0, D0 = get_inputs_and_disturbances(cgina_val, Fgina_val, F_val, Q_val, N)

    F0 = U0[1,:]
    lb = fill(0.22, N) # Minimum F
    ub = fill(0.4, N) # Maximum F

    # Define ODE problem
    prob = ODEProblem(system_equations!, x0, (t0, tf), p)

    # Define objective function for specific scenario
    function obj_scenario(F, p) 
        return loss_cap_eff_ref(hcat(F, U0[2,:])', D0, prob, t_vec, cap_eff_ref, lambda_dU)
    end

    # Define constraints
    cons = nothing
    lb_con = nothing
    ub_con = nothing

    # Collect all scenario-specific data into a tuple to return
    scenario_data = (obj_scenario = obj_scenario, U0 = U0, D0 = D0, cap_eff_ref = cap_eff_ref, lb = lb, ub = ub, prob = prob, t_vec = t_vec, lambda_dU = lambda_dU, x0 = x0)

    return obj_scenario, F0, lb, ub, cons, lb_con, ub_con, scenario_data
end


function define_scenario01B(t0, tf, p, N) 
    ### Scenario 2 ###

    t_vec = range(t0, tf, length=N+1)

    # Define input and disturbance sequence for scenario
    cgina_val = [4.4] 
    Fgina_val = [180]
    F_val = [0.26]
    Q_val = [30,25]
    cap_eff_ref_val = [90]
    cap_eff_ref = repeat(cap_eff_ref_val, inner=Int(N/length(cap_eff_ref_val)))
    # MPC tuning parameters
    lambda_dU = 1e1

    # Initial state
    x0 = [2.63; 0.3; 0.3]

    U0, D0 = get_inputs_and_disturbances(cgina_val, Fgina_val, F_val, Q_val, N)

    F0 = U0[1,:]
    lb = fill(0.22, N) # Minimum F
    ub = fill(0.4, N) # Maximum F

    # Define ODE problem
    prob = ODEProblem(system_equations!, x0, (t0, tf), p)

    # Define objective function for specific scenario
    function obj_scenario(F, p) 
        return loss_cap_eff_ref(hcat(F, U0[2,:])', D0, prob, t_vec, cap_eff_ref, lambda_dU)
    end

    # Define constraints
    cons = nothing
    lb_con = nothing
    ub_con = nothing

    # Collect all scenario-specific data into a tuple to return
    scenario_data = (obj_scenario = obj_scenario, U0 = U0, D0 = D0, cap_eff_ref = cap_eff_ref, lb = lb, ub = ub, prob = prob, t_vec = t_vec, lambda_dU = lambda_dU, x0 = x0)

    return obj_scenario, F0, lb, ub, cons, lb_con, ub_con, scenario_data
end


function define_scenario01C(t0, tf, p, N) 
    ### Scenario 3 ###

    t_vec = range(t0, tf, length=N+1)

    # Define input and disturbance sequence for scenario
    cgina_val = rand(N) * 0.6 .+ 4.1 # Random values between 4.3 and 4.5
    Fgina_val = [170]
    F_val = [0.26]
    Q_val = [27]
    cap_eff_ref_val = [90, 84, 94, 88] # Varying reference trajectory
    cap_eff_ref = repeat(cap_eff_ref_val, inner=Int(N/length(cap_eff_ref_val)))
    # MPC tuning parameters
    lambda_dU = 1e4

    # Initial state
    x0 = [2.63; 0.3; 0.3]

    U0, D0 = get_inputs_and_disturbances(cgina_val, Fgina_val, F_val, Q_val, N)

    F0 = U0[1,:]
    lb = fill(0.22, N) # Minimum F
    ub = fill(0.4, N) # Maximum F

    # Define ODE problem
    prob = ODEProblem(system_equations!, x0, (t0, tf), p)

    # Define objective function for specific scenario
    function obj_scenario(F, p) 
        return loss_cap_eff_ref(hcat(F, U0[2,:])', D0, prob, t_vec, cap_eff_ref, lambda_dU)
    end

    # Define constraints
    cons = nothing
    lb_con = nothing
    ub_con = nothing

    # Collect all scenario-specific data into a tuple to return
    scenario_data = (obj_scenario = obj_scenario, U0 = U0, D0 = D0, cap_eff_ref = cap_eff_ref, lb = lb, ub = ub, prob = prob, t_vec = t_vec, lambda_dU = lambda_dU, x0 = x0)

    return obj_scenario, F0, lb, ub, cons, lb_con, ub_con, scenario_data
end


# EMPC scenarios
function define_scenario02A(t0, tf, p, N) 
    ### Scenario 2A ###
    # Economic MPC

    t_vec = range(t0, tf, length=N+1)

    # Define input and disturbance sequence for scenario
    cgina_val = [4.4]
    Fgina_val = [170]
    F_val = [0.3]
    Q_val = [27]

    # Minimum average capture efficiency
    eta_min = 90.0
    
    # MPC tuning parameters
    lambda = [1e4, 1e4, 1e4]

    # Initial state
    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # ca, cd, ct, Ma, Mina

    U0_bar, D0 = get_inputs_and_disturbances(cgina_val, Fgina_val, F_val, Q_val, N)

    U0 = vcat(vec(U0_bar'), [0]) # Append initial slack variable value to input sequence

    # El. prices and CO2 prices
    p_el = repeat([50, 100], inner=Int(N/2))
    p_CO2 = repeat([10], inner=N)
    D0_w_prices = vcat(D0, p_el', p_CO2')

    # Lower and upper bounds for optimization variables
    Fmin = 0.22
    Fmax = 0.4
    Qmin = 22
    Qmax = 35
    smin = 0
    smax = Inf
    lb = vcat(fill(Fmin, N), fill(Qmin, N), smin)
    ub = vcat(fill(Fmax, N), fill(Qmax, N), smax)

    # Define ODE problem
    prob = ODEProblem(system_equations_w_cum_CO2!, x0, (t0, tf), p)

    # Define objective function for specific scenario
    function obj_scenario(U, p)
        return loss_EMPC(U, D0_w_prices, prob, t_vec, lambda, eta_min)  
        
    end

    # Define constraints
    function cons(res, U, p)
        # Extract inputs
        F = U[1:N]
        Q = U[N .+ (1:N)]
        s = U[2N + 1]

        Usim = hcat(F, Q)' # Combine F and Q into input matrix for simulation

        X, Z = simulate_control_intervals(Usim, D0_w_prices, prob, t_vec)

        Ma_tf = X[4,end] # Cumulative captured CO2 at final time
        Mina_tf = X[5,end] # Cumulative inlet CO2 at final time
        
        res .= Ma_tf - Mina_tf * eta_min / 100.0 + s

    end

    lb_con = [0]
    ub_con = [Inf]

    # Collect all scenario-specific data into a tuple to return
    scenario_data = (obj_scenario = obj_scenario, U0 = U0, D0 = D0_w_prices, eta_min = eta_min, lb = lb, ub = ub, prob = prob, t_vec = t_vec, lambda = lambda, x0 = x0)

    return obj_scenario, U0, lb, ub, cons, lb_con, ub_con, scenario_data
end