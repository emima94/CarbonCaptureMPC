function define_scenario01A()

    t0 = 0
    tf = 6
    dt = 0.25
    N = Int((tf - t0) / dt)  # number of intervals
    t = range(t0, tf, length=N+1)

    ODE_solver = "euler"
    #ODE_solver = "rk4"

    Nx_per_interval = 10
    Nx = N * Nx_per_interval + 1   # total number of state nodes
    dt_sim = dt / Nx_per_interval
    t_sim = range(t0, tf, length=Nx)

    # What control interval does each state belong to?
    interval_of = repeat(1:N, inner=Nx_per_interval)
    push!(interval_of, N)  # last state belongs to last interval

    # Define scenario
    cgina_val = [4.4]
    Fgina_val = [160, 180]
    cgina = repeat(cgina_val, inner = Int(N/length(cgina_val))) # Disturbance sequence (cgina)
    Fgina = repeat(Fgina_val, inner = Int(N/length(Fgina_val))) # Disturbance sequence (Fgina)
    D = hcat(cgina, Fgina)'
    Q = repeat([28.0], inner=N) # Disturbance sequence (Q)

    # Previous control input (will be updated in closed-loop)
    F_prev = 0.3

    # Minimum capture efficiency
    cap_eff_ref_val = [90]
    cap_eff_ref = repeat(cap_eff_ref_val, inner = Int((Nx-1)/length(cap_eff_ref_val)))

    # MPC tuning Parameters
    lambda_F = 1e4
    lambda = lambda_F

    # ODE rhs:
    #ffun = system_equations_w_cum_CO2
    ffun = system_equations_sys_id_w_cum_CO2
    # Output equations
    #hfun = output_equations
    hfun = output_equations_sys_id_detailed

    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state

    # Bounds 
    F_min, F_max = 0.22, 0.4
    lb = F_min
    ub = F_max

    # Scenario name 
    name = "scenario01A"

    p01A = (N = N, Nx = Nx, interval_of = interval_of, D = D, Q = Q, dt = dt, dt_sim = dt_sim, t = t, t_sim = t_sim, ffun = ffun, hfun = hfun, F_prev = F_prev, cap_eff_ref = cap_eff_ref, lambda = lambda, x0 = x0, lb = lb, ub = ub, name = name, ODE_solver = ODE_solver)
    return p01A

end

function define_scenario01B()

    t0 = 0
    tf = 6
    dt = 0.25
    N = Int((tf - t0) / dt)  # number of intervals
    t = range(t0, tf, length=N+1)

    # ODE_solver
    ODE_solver = "rk4"

    Nx_per_interval = 5
    Nx = N * Nx_per_interval + 1   # total number of state nodes
    dt_sim = dt / Nx_per_interval
    t_sim = range(t0, tf, length=Nx)

    # What control interval does each state belong to?
    interval_of = repeat(1:N, inner=Nx_per_interval)
    push!(interval_of, N)  # last state belongs to last interval

    # Define scenario
    cgina_val = [4.4]
    Fgina_val = [180]
    cgina = repeat(cgina_val, inner = Int(N/length(cgina_val))) # Disturbance sequence (cgina)
    Fgina = repeat(Fgina_val, inner = Int(N/length(Fgina_val))) # Disturbance sequence (Fgina)
    D = hcat(cgina, Fgina)'
    Q = repeat([30, 25], inner=Int(N/2)) # Disturbance sequence (Q)

    # Previous control input (will be updated in closed-loop)
    F_prev = 0.3

    # Minimum capture efficiency
    cap_eff_ref_val = [90.0]
    cap_eff_ref = repeat(cap_eff_ref_val, inner = Int((Nx-1)/length(cap_eff_ref_val)))

    # MPC tuning Parameters
    #lambda_F = 1e3
    lambda_F = 1e4
    lambda = lambda_F

    # ODE rhs:
    ffun = system_equations_w_cum_CO2
    # Output equations
    hfun = output_equations

    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state

    # Bounds 
    F_min, F_max = 0.22, 0.4
    lb = F_min
    ub = F_max

    # Scenario name 
    name = "scenario01B"

    p01B = (N = N, Nx = Nx, interval_of = interval_of, D = D, Q = Q, dt = dt, dt_sim = dt_sim, t = t, t_sim = t_sim, ffun = ffun, hfun = hfun, F_prev = F_prev, cap_eff_ref = cap_eff_ref, lambda = lambda, x0 = x0, lb = lb, ub = ub, name = name, ODE_solver = ODE_solver)
    return p01B

end


function define_scenario01C()

    t0 = 0
    tf = 6
    dt = 0.25
    N = Int((tf - t0) / dt)  # number of intervals
    t = range(t0, tf, length=N+1)

    # ODE_solver
    ODE_solver = "euler"

    Nx_per_interval = 50
    Nx = N * Nx_per_interval + 1   # total number of state nodes
    dt_sim = dt / Nx_per_interval
    t_sim = range(t0, tf, length=Nx)

    # What control interval does each state belong to?
    interval_of = repeat(1:N, inner=Nx_per_interval)
    push!(interval_of, N)  # last state belongs to last interval

    # Define scenario
    cgina_val = rand(N) * 0.5 .+ 4.2 # Disturbance sequence (cgina) with some variability
    Fgina_val = [180]
    cgina = repeat(cgina_val, inner = Int(N/length(cgina_val))) # Disturbance sequence (cgina)
    Fgina = repeat(Fgina_val, inner = Int(N/length(Fgina_val))) # Disturbance sequence (Fgina)
    D = hcat(cgina, Fgina)'
    Q = repeat([27], inner=N) # Disturbance sequence (Q)

    # Previous control input (will be updated in closed-loop)
    F_prev = 0.3

    # Minimum capture efficiency
    cap_eff_ref_val = [85.0, 95.0, 82.0, 88.0]
    cap_eff_ref = repeat(cap_eff_ref_val, inner = Int((Nx-1)/length(cap_eff_ref_val)))

    # MPC tuning Parameters
    #lambda_F = 500
    lambda_F = 1e4
    lambda = lambda_F

    # ODE rhs:
    ffun = system_equations_w_cum_CO2
    # Output equations
    hfun = output_equations

    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state

    # Bounds 
    F_min, F_max = 0.22, 0.4
    lb = F_min
    ub = F_max

    # Scenario name 
    name = "scenario01C"

    p01C = (N = N, Nx = Nx, interval_of = interval_of, D = D, Q = Q, dt = dt, dt_sim = dt_sim, t = t, t_sim = t_sim, ffun = ffun, hfun = hfun, F_prev = F_prev, cap_eff_ref = cap_eff_ref, lambda = lambda, x0 = x0, lb = lb, ub = ub, name = name, ODE_solver = ODE_solver)
    return p01C

end


function define_scenario02A()

    t0 = 0
    tf = 12
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
    eta_min = 0.80

    # Electricity and CO2 prices
    p_el_val = [0.2, 0.2, 1.0, 1.0, 0.2, 0.2] # EUR/kWh
    p_el = repeat(p_el_val, inner=Int(N/length(p_el_val)))
    p_CO2 = 0.07 * 10  #EUR/kg

    # MPC tuning Parameters
    lambda_F = 1e3
    lambda_Q = 1e-1
    lambda_s = 1e6
    lambda = [lambda_F, lambda_Q, lambda_s]

    # ODE rhs:
    ffun = system_equations_w_cum_CO2
    # Output equations
    hfun = output_equations

    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state

    # Bounds 
    F_min, F_max = 0.22, 0.32
    Q_min, Q_max = 22.0, 35.0
    s_min = 0.0

    lb = (F_min = F_min, Q_min = Q_min, s_min = s_min)
    ub = (F_max = F_max, Q_max = Q_max, s_max = Inf)

    p02A = (p_el = p_el, p_CO2 = p_CO2, lambda = lambda, eta_min = eta_min, x0 = x0, lb = lb, ub = ub, N = N, Nx = Nx, interval_of = interval_of, D = D, dt = dt, dt_sim = dt_sim, t = t, t_sim = t_sim, ffun = ffun, hfun = hfun, F_prev = F_prev, Q_prev = Q_prev)

    return p02A

end


function define_scenario_sim()

    t0 = 0.0
    tf = 24.0
    dt = 0.25
    N = Int((tf - t0) / dt)  # number of intervals
    t = range(t0, tf, length=N+1)

    Nx_per_interval = 50
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

    #F_val = [0.25, 0.35, 0.3, 0.4, 0.22, 0.30]
    #Q_val = [25, 30, 28, 35, 22, 27]
    F_val = [0.4]
    Q_val = [26.0]
    U_sim = hcat(repeat(F_val, inner=Int(N/length(F_val))), repeat(Q_val, inner=Int(N/length(Q_val))))'

    # Previous control input (will be updated in closed-loop)
    F_prev = 0.3
    Q_prev = 28.0

    # Minimum capture efficiency
    eta_min = 0.80

    # Electricity and CO2 prices
    p_el_val = [0.2, 0.2, 1.0, 1.0, 0.2, 0.2] # EUR/kWh
    p_el = repeat(p_el_val, inner=Int(N/length(p_el_val)))
    p_CO2 = 0.07  #EUR/kg

    # MPC tuning Parameters
    lambda_F = 1e2
    lambda_Q = 1e-1
    lambda_s = 1e6
    lambda = [lambda_F, lambda_Q, lambda_s]

    # ODE rhs:
    ffun = system_equations_w_cum_CO2
    # Output equations
    hfun = output_equations

    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state

    # Bounds 
    F_min, F_max = 0.22, 0.4
    Q_min, Q_max = 22.0, 35.0
    s_min = 0.0

    lb = (F_min = F_min, Q_min = Q_min, s_min = s_min)
    ub = (F_max = F_max, Q_max = Q_max, s_max = Inf)

    p02A = (p_el = p_el, p_CO2 = p_CO2, lambda = lambda, eta_min = eta_min, x0 = x0, lb = lb, ub = ub, N = N, Nx = Nx, interval_of = interval_of, D = D, dt = dt, dt_sim = dt_sim, t = t, t_sim = t_sim, ffun = ffun, hfun = hfun, F_prev = F_prev, Q_prev = Q_prev, U_sim = U_sim)

    return p02A

end



function define_scenario03()

    # Closed loop scenario
    t0 = 0
    tf = 12
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
    p_el_val = [0.1, 0.12, 0.12, 0.24, 0.26, 0.31, 0.25, 0.15, 0.02, 0.02, 0.05, 0.03, 0.16, 0.44, 0.65, 0.56, 0.65, 0.7, 0.10, 0.04, 0.02, 0.10, 0.09, 0.08] # EUR/kWh
    p_el = repeat(p_el_val, inner=Int(N/length(p_el_val)))
    p_CO2 = 0.07 * 10  #EUR/kg

    # MPC tuning Parameters
    lambda_F = 1e3
    lambda_Q = 1e-1
    lambda_s = 1e6
    lambda = [lambda_F, lambda_Q, lambda_s]

    # ODE rhs:
    ffun = system_equations_w_cum_CO2
    # Output equations
    hfun = output_equations

    x0 = [2.63; 0.3; 0.3; 0.0; 0.0] # Initial state

    # Bounds 
    F_min, F_max = 0.22, 0.32
    Q_min, Q_max = 22.0, 35.0
    s_min = 0.0

    lb = (F_min = F_min, Q_min = Q_min, s_min = s_min)
    ub = (F_max = F_max, Q_max = Q_max, s_max = Inf)

    p02A = (p_el = p_el, p_CO2 = p_CO2, lambda = lambda, eta_min = eta_min, x0 = x0, lb = lb, ub = ub, N = N, Nx = Nx, interval_of = interval_of, D = D, dt = dt, dt_sim = dt_sim, t = t, t_sim = t_sim, ffun = ffun, hfun = hfun, F_prev = F_prev, Q_prev = Q_prev)

    return p02A

end