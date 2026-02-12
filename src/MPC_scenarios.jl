
function define_scenario01(t0, tf, p, N) 
    ### Scenario 1 ###

    t_vec = range(t0, tf, length=N+1)

    # Define input and disturbance sequence for scenario
    cgina_val = [4.4]
    Fgina_val = [160, 165, 165, 150, 150, 170, 170, 160]
    F_val = [0.3]
    Q_val = [36,24]
    cap_eff_ref_val = [85, 95, 80, 90]
    cap_eff_ref = repeat(cap_eff_ref_val, inner=Int(N/length(cap_eff_ref_val)))
    # MPC tuning parameters
    lambda_dU = 1e2

    # Initial state
    x0 = [2.63; 0.3; 0.3]

    U0, D0 = get_inputs_and_disturbances(cgina_val, Fgina_val, F_val, Q_val, N)

    F0 = U0[1,:]
    lb = fill(0.1, N) # Minimum F
    ub = fill(0.4, N) # Maximum F

    # Define ODE problem
    prob = ODEProblem(system_equations!, x0, (t0, tf), p)

    # Define objective function for specific scenario
    function obj_scenario01(F, p) 
        return loss_cap_eff_ref(hcat(F, U0[2,:])', D0, prob, t_vec, cap_eff_ref, lambda_dU)
    end

    # Collect all scenario-specific data into a tuple to return
    scenario_data = (obj_scenario01 = obj_scenario01, U0 = U0, D0 = D0, cap_eff_ref = cap_eff_ref, lb = lb, ub = ub, prob = prob, t_vec = t_vec, lambda_dU = lambda_dU, x0 = x0)

    return obj_scenario01, F0, lb, ub, scenario_data
end