
# Define system equations
function system_equations!(dx, x, p, t)
    
    # Extract states
    ca, cd, ct = x

    # Extract from p
    u, d, theta = p

    # Parameters #
    NA_params, ND_params, s, V, out_allocate = theta
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Unit conversion
    mol2kmol = 1e-3

    # Evaluate Na
    input_Na = [ca]
    Na = evaluate_NN(NA_params, input_Na)[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = evaluate_NN(ND_params, input_Nd)[1] # Output is scalar

    # ODEs:
    dx[1] = (F * (ct - ca) + Na * mol2kmol) / Va
    dx[2] = (F * (ca - cd) - Nd * mol2kmol) / Vd
    dx[3] = F * (cd - ct) / Vt
    
    return nothing

end

function system_equations_w_cum_CO2!(dx, x, p, t)

    # Extract states
    ca, cd, ct, Ma, Mina = x

    # Extract from p
    u, d, theta = p

    # Parameters #
    NA_params, ND_params, s, V, out_allocate = theta
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Conversion Constants
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Inputs #
    F, Q = u

    # Disturbances #
    cgina, Fgina = d

    # Evaluate Na
    input_Na = [ca]
    Na = evaluate_NN(NA_params, input_Na)[1] # Output is scalar

    # Evaluate Nd
    input_Nd = [ca, cd, F, Q]
    Nd = evaluate_NN(ND_params, input_Nd)[1] # Output is scalar

    # ODEs:
    dx[1] = (F * (ct - ca) + Na * mol2kmol) / Va
    dx[2] = (F * (ca - cd) - Nd * mol2kmol) / Vd
    dx[3] = F * (cd - ct) / Vt
    dx[4] = Na * M_CO2 * g2kg # Cumulative CO2 captured
    dx[5] = Fgina * cgina * M_CO2 * g2kg # Cumulative CO2 in
    
    return nothing

end

function output_equations(x, p, t)
    # Extract states
    ca, _, _ = x

    # Extract from p
    u, d, theta = p

    # Parameters #
    NA_params, ND_params, s, V, out_allocate = theta
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Inputs #
    F, Q = u
    
    # Disturbances #
    cgina, Fgina = d

    # Constants 
    R = 8.314 # J/(mol*K)
    Ta = 313.15 # K
    Pa = 101325 # Pa
    M_CO2 = 44.01 # g/mol
    g2kg = 1e-3 # g to kg
    mol2kmol = 1e-3

    # Evaluate Na
    input_Na = [ca]
    Na = evaluate_NN(NA_params, input_Na)[1] # Output is scalar

    # Mass balance, CO2, gas phase, absorber, steady state:
    cga = cgina - Na / Fgina
    
    # Total mol concentration in gas phase (assuming ideal gas behavior and constant T and P)
    cga_tot = Pa / (R * Ta)

    # Inlet CO2 concentration [vol%]:
    ygina = cgina / cga_tot * 100

    # Output CO2 concentration [vol%]:
    yga = cga / cga_tot * 100

    # Capture efficiency:
    cap_eff = (ygina - yga) / ygina * 100

    # Capture rate [kg/h]
    ma = Na * M_CO2 * g2kg

    z = [yga, ygina, cap_eff, ma]

    return z
end


