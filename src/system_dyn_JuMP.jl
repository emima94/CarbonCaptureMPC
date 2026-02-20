function system_equations_w_cum_CO2(x, u, d, p)

    # Extract states
    ca, cd, ct, _, _ = x

    # Parameters #
    NA_params, ND_params, s, V = p
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
    dx = [(F * (ct - ca) + Na * mol2kmol) / Va,
        (F * (ca - cd) - Nd * mol2kmol) / Vd,
        F * (cd - ct) / Vt,
        Na * M_CO2 * g2kg, # Cumulative CO2 captured
        Fgina * cgina * M_CO2 * g2kg # Cumulative CO2 in
    ] 
    
    return dx

end

function system_equations(x, u, d, p)

    # Extract states
    ca, cd, ct = x

    # Parameters #
    NA_params, ND_params, s, V = p
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
    dx = [(F * (ct - ca) + Na * mol2kmol) / Va,
        (F * (ca - cd) - Nd * mol2kmol) / Vd,
        F * (cd - ct) / Vt
    ] 
    
    return dx

end


function output_equations(x, u, d, p)
    # Extract states
    ca = x[1]
    cd = x[2]

    # Parameters #
    NA_params, ND_params, s, V = p
    sa, sd, st, s_yga = s
    Va, Vd, Vt = V

    # Conversion Constants
    R = 8.314 # J/(mol*K)
    Ta = 313.15 # K
    Pa = 101325 # Pa
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

    z = [yga, ygina, cap_eff, ma, Na, Nd]

    return z

end