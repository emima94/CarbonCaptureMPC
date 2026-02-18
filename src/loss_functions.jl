

function loss_cap_eff_ref(U, D, prob, t_vec, cap_eff_ref, lambda_dU)

    """
    Reference control of capture efficiency.
    MV: F
    cost: (cap_eff - cap_eff_ref)^2 + lambda_dU * sum(diff(F).^2)

    """
    _, Z = simulate_control_intervals(U, D, prob, t_vec)
    
    cap_eff = Z[3,:] 

    # Compute loss (e.g., squared error from target trajectory)
    loss_ref = sum((cap_eff .- cap_eff_ref).^2)
    loss_dU = lambda_dU * sum(diff(U[1,:]).^2) # Penalize large changes in F
    loss = loss_ref + loss_dU
    return loss

end

function loss_EMPC(U, D, prob, t_vec, lambda, eta_min)
    """
    Economic MPC objective function.
    MV: F, Q, s

    """
    N = length(t_vec) - 1
    dt = diff(t_vec)

    # Extract inputs
    F = U[1:N]
    Q = U[N .+ (1:N)]
    s = U[2N + 1]

    # Extract disturbances
    cgina = D[1, :]
    Fgina = D[2, :]
    p_el = D[3, :]
    p_CO2 = D[4, :]

    # Inlet CO2 mole flow rate [mol/h]
    m_ina = cgina .* Fgina

    # Extract tuning parameters
    lambda_F = lambda[1]
    lambda_Q = lambda[2]
    lambda_s = lambda[3]

    U_sim = hcat(F, Q)'
    D_sim = hcat(cgina, Fgina)'
    
    _, Z = simulate_control_intervals(U_sim, D_sim, prob, t_vec)

    # Extract output
    ma = Z[4,:]     # captured CO2 [kg/h]
    
    cost = sum(p_el .* Q .* dt) - sum(p_CO2 .* ma .* dt) + lambda_F * sum(diff(F).^2) + lambda_Q * sum(diff(Q).^2) + lambda_s * s^2

    return cost
end

