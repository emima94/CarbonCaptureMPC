
# Define system equations
function system_equations!(dx, x, p, t)
    
    # Extract states
    ca, cd, ct = x

    # Extract inputs, disturbances, and parameters
    u, d, theta = p
    F, Q = u
    cgina, Fgina = d
    NA_params, ND_params, s, V = theta
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Evaluate Na
    input_Na = [ca]
    Na = evaluate_NN(NA_params, input_Na)[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = evaluate_NN(ND_params, input_Nd)[1] # Output is scalar

    # ODEs:
    dx[1] = (F * (ct - ca) + Na) / Va
    dx[2] = (F * (ca - cd) - Nd) / Vd
    dx[3] = F * (cd - ct) / Vt
    
    return nothing

end