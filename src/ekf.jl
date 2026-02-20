function ekf_update(y_meas, x_prior, y_prior, P_prior, R_prior, H_k, p, t_cur, obs_eq, S)
        
    n_states = length(x_prior)

    # Innovation (residual)
    e_k = y_meas - y_prior

    # Kalman gain
    K_k = P_prior * H_k' * (R_prior \ I)  # equivalently: P_prior[k]*H_k'*inv(R_prior[k])

    # Updated state estimate
    x_post = x_prior + K_k * e_k

    # Updated state covariance
    I_n = Matrix{Float64}(I, n_states, n_states)
    P_post = (I_n - K_k * H_k) * P_prior

    # Updated observation (posterior)
    y_post = obs_eq(x_post, p, t_cur)

    R_post = H_k * P_post * H_k' + S

    return x_post, y_post, P_post, R_post

end