

function plot_EMPC(model, p_EMPC, p, U_sim = nothing)
    # Extract scenario parameters
    p_el = p_EMPC.p_el
    p_CO2 = p_EMPC.p_CO2
    lambda = p_EMPC.lambda
    eta_min = p_EMPC.eta_min
    x0 = p_EMPC.x0
    lb = p_EMPC.lb
    ub = p_EMPC.ub
    N = p_EMPC.N
    Nx = p_EMPC.Nx
    Nx_per_interval = Int((Nx-1) / N)
    interval_of = p_EMPC.interval_of
    D = p_EMPC.D
    dt_sim = p_EMPC.dt_sim
    t = p_EMPC.t
    t_sim = p_EMPC.t_sim
    ffun = p_EMPC.ffun
    hfun = p_EMPC.hfun
    F_prev = p_EMPC.F_prev
    Q_prev = p_EMPC.Q_prev

    lambda_F, lambda_Q, lambda_s = lambda

    # Bounds
    F_min, Q_min, s_min = lb
    F_max, Q_max, s_max = ub


    # Extract values
    # x_opt = value.(model[:x])
    # F_opt = value.(model[:F])
    # Q_opt = value.(model[:Q])
    # s_opt = value(model[:s])

    if U_sim === nothing
         F_sim = value.(model[:F])
         Q_sim = value.(model[:Q])
    else
        F_sim = U_sim[1,:]
        Q_sim = U_sim[2,:]
    end

    U_sim = repeat(hcat(F_sim, Q_sim)', inner=(1,Nx_per_interval))
    Dsim = repeat(D, inner=(1,Nx_per_interval))

    X, Z = euler_simulation(ffun, hfun, x0, U_sim, Dsim, p, t_sim, dt_sim)


    # Capture rate
    ma = Z[4,:]     # captured CO2 [kg/h]

    # Accumulated CO2 captured
    Ma = X[4,:] # Cumulative CO2 captured [kg]

    # Compute cost and revenues
    cost = p_el .* Q_sim # Cost of electricity [EUR/h]
    revenue = p_CO2 * ma # [EUR/h] Revenue from captured CO2
    profit = revenue .- repeat(cost, inner=Nx_per_interval) # Total profit over time [EUR/h]

    # Compute SRD
    s2h = 1/3600 # convert from per second to per hour
    SRD = (repeat(Q_sim, inner=Nx_per_interval) / s2h) ./ ma * 1e-3 # kJ/h per kg/h of CO2 captured

    colors = palette(:auto, 10)

    # Plotting
    plt_cap_eff = plot(t_sim[2:end], Z[3,:], seriestype=:steppre, label=L"\eta", title="Output", xlabel="", ylabel="Cap. eff. [%]", ylims=(50,110))


    plt2 = plot(t_sim[2:end], Z[1:2,:]', labels=["yga" "ygina"], title="Output Trajectories with Optimized F", xlabel ="", ylabel="Concentration (vol%)", ylims = (-2,14))
    plt3 = plot(t_sim[2:end], Z[4,:], labels=["Na"],  title = "Capture rate", xlabel ="", ylabel="Capture Rate (kg/hr)", ylims=(0,50))

    plt_x = plot(t_sim, X[1:3,:]', labels=[L"c_a" L"c_d" L"c_t"], title="States", xlabel="", ylabel="CO2 liq. conc. [kmol/m3]", ylims = (-1,5))
    lambda_F_str = @sprintf("%.0e", lambda_F)
    lambda_Q_str = @sprintf("%.0e", lambda_Q)

    plt_F = plot(t, vcat(F_sim[1], F_sim), seriestype = :steppre, title=L"\lambda_F = %$lambda_F_str", xlabel ="", ylabel="Flow Rate", ylims=(0.2,0.35), label=L"F")
    plot!(plt_F, t, fill(F_min, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
    plot!(plt_F, t, fill(F_max, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")

    plt_Q = plot(t, vcat(Q_sim[1], Q_sim), seriestype = :steppre, label=L"Q", title=L"\lambda_Q = %$lambda_Q_str", xlabel="", ylabel="Reboiler duty [kW]", ylims = (20,37))
    plot!(plt_Q, t, fill(Q_min, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
    plot!(plt_Q, t, fill(Q_max, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")

    plt_cgina = plot(t_sim, vcat(Dsim[1,1], Dsim[1,:]), seriestype = :steppre, label=L"c_{in,a}^g", title="Disturbances", xlabel="", ylabel="Flue gas conc. [mol/m3]", ylims=(4,5))
    plt_Fgina = plot(t_sim, vcat(Dsim[2,1], Dsim[2,:]), seriestype = :steppre, label=L"F_{in,a}^g", title="", xlabel="", ylabel="Flue gas flow rate [m3/h]", ylims=(140,220))

    plt_p_el = plot(t, vcat(p_el[1], p_el) * 1000, seriestype=:steppre, label=L"p_{el}", title="Electricity Price", xlabel="", ylabel="Price [EUR/MWh]", ylims=(0,1200))
    plt_profit = plot(t_sim, vcat(profit[1], profit), seriestype=:steppre, label="", title="Profit", xlabel ="", ylabel="Profit [EUR/h]")

    plt_cap_cum = plot(t_sim, Ma, seriestype=:steppre, label=L"M_a", title="Cumulative CO2 Captured", xlabel ="", ylabel="Mass of CO2 Captured [kg]", ylims=(0,600))
    hline!([eta_min * X[5,end]], linestyle=:dash, color=:red, label="")
    annotate!(t_sim[end]*0.5, eta_min * X[5,end] * 1.05, text("Req.: $(eta_min * 100) %", font("Times", :black, 8)))
    hline!([X[5,end]], linestyle=:dash, color=:blue, label="")
    annotate!(t_sim[end]*0.5, X[5,end] * 1.05, text("Max.", font("Times", :black, 8)))

    plt_SRD = plot(t_sim, vcat(SRD[1], SRD), seriestype=:steppre, label=L"SRD", title="Specific Reboiler Duty", xlabel ="", ylabel="SRD [MJ/kg]", ylims=(0,5))

    scale = 1200.0
    H2W_ratio = 0.4

    # Write title as annotation in the first plot
    total_profit = sum(profit) * dt_sim
    title_str = "Total profit: $(round(total_profit, digits=1)) EUR"
    #annotate!(plt_cap_eff, 0.5 * t_sim[end], 120, text(title_str, font("Times", :black, 10, :bold)))

    plt_out = plot(plt_cap_eff, plt_cap_cum, plt_SRD, plt_cgina, plt_F, 
                    plt_x, plt_profit, plt_p_el, plt_Fgina, plt_Q,
        layout = (2,5), link=:x, 
        size = (scale, scale * H2W_ratio),
        left_margin = 6mm,
        plot_title = title_str
    )

    return plt_out, total_profit
end


function plot_ref_MPC(model, p_ref, p)
    # Extract scenario parameters for reference control scenario
    ffun = p_ref.ffun
    hfun = p_ref.hfun
    x0 = p_ref.x0
    D = p_ref.D
    t_sim = p_ref.t_sim
    t = p_ref.t
    dt_sim = p_ref.dt_sim
    F_prev = p_ref.F_prev
    cap_eff_ref = p_ref.cap_eff_ref
    lambda_F = p_ref.lambda
    Q = p_ref.Q
    N = p_ref.N
    Nx = p_ref.Nx
    interval_of = p_ref.interval_of
    Nx_per_interval = Int((Nx-1) / N)
    F_min = p_ref.lb
    F_max = p_ref.ub

    # Extract values
    x_opt = value.(model[:x])
    F_opt = value.(model[:F])

    U_sim = repeat(hcat(F_opt, Q)', inner=(1,Nx_per_interval))
    Dsim = repeat(D, inner=(1,Nx_per_interval))

    X, Z = euler_simulation(ffun, hfun, x0, U_sim, Dsim, p, t_sim, dt_sim)

    # Capture rate
    ma = Z[4,:]     # captured CO2 [kg/h]

    # Accumulated CO2 captured
    Ma = X[4,:] # Cumulative CO2 captured [kg]

    # Compute SRD
    s2h = 1/3600 # convert from per second to per hour
    SRD = (repeat(Q, inner=Nx_per_interval) / s2h) ./ ma * 1e-3 # kJ/h per kg/h of CO2 captured

    colors = palette(:auto)

    # Plotting
    plt_cap_eff = plot(t_sim[2:end], cap_eff_ref, seriestype=:steppre, color = colors[2], linewidth = 2, linestyle=:dot, label=L"\eta_{ref}", title="Output", xlabel="", ylabel="Cap. eff. [%]", ylims=(70,110))
    plot!(plt_cap_eff, t_sim[2:end], Z[3,:], label=L"\eta", color = colors[1])

    plt_x = plot(t_sim, X[1:3,:]', labels=[L"c_a" L"c_d" L"c_t"], title="States", xlabel="", ylabel="CO2 liq. conc. [kmol/m3]", ylims = (-1,5))
    lambda_F_str = @sprintf("%.0e", lambda_F)
    plt_F = plot(t, vcat(F_opt[1], F_opt), seriestype = :steppre, title=L"\lambda_F = %$lambda_F_str", xlabel ="", ylabel="Flow Rate [m3/h]", ylims=(0.2,0.45), label=L"F")
    plot!(plt_F, t, fill(F_min, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
    plot!(plt_F, t, fill(F_max, N+1), seriestype=:steppre, linestyle=:dash, color=:red, label = "")

    plt_Q = plot(t, vcat(Q[1], Q), seriestype = :steppre, label=L"Q", title="", xlabel="", ylabel="Reboiler duty [kW]", ylims = (20,37))
    
    plt_cgina = plot(t_sim, vcat(Dsim[1,1], Dsim[1,:]), seriestype = :steppre, label=L"c_{in,a}^g", title="Disturbances", xlabel="", ylabel="Flue gas conc. [mol/m3]", ylims=(4,5))
    plt_Fgina = plot(t_sim, vcat(Dsim[2,1], Dsim[2,:]), seriestype = :steppre, label=L"F_{in,a}^g", title="", xlabel="", ylabel="Flue gas flow rate [m3/h]", ylims=(140,220))

    scale = 800.0
    H2W_ratio = 0.6

    # Add x-axis to specified plots
    plt_x_axis_names = [plt_x, plt_Fgina, plt_Q]
    for plt in plt_x_axis_names
        plot!(plt, xlabel="Time [h]")
    end

    plt_out = plot(plt_cap_eff, plt_cgina, plt_F, 
                    plt_x, plt_Fgina, plt_Q,
        layout = (2,3), link=:x, 
        size = (scale, scale * H2W_ratio),
        left_margin = 6mm
    )

    return plt_out
end