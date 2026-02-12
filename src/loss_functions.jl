

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