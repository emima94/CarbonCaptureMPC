## Plots 

# Plot results
function plot_results(Usim, Dsim, Xsim, Zsim, t_fine_sim, scen_data)
    plt_cap_eff = plot(t_vec, vcat(scen_data.cap_eff_ref[1], scen_data.cap_eff_ref), seriestype=:steppre, linestyle=:dot, label=L"y^g_{a,ref}", title="Output", xlabel="", ylabel="Cap. eff. [%]", ylims=(75,110))
    plot!(plt_cap_eff, t_fine_sim[2:end], Zsim[3,:], label=L"y^g_a")

    plt2 = plot(t_fine_sim[2:end], Zsim[1:2,:]', labels=["yga" "ygina"], title="Output Trajectories with Optimized F", xlabel="Time (hr)", ylabel="Concentration (vol%)", ylims = (-2,14))
    plt3 = plot(t_fine_sim[2:end], Zsim[4,:], labels=["Na"],  title = "Capture rate", xlabel="Time (hr)", ylabel="Capture Rate (mol/hr)", ylims=(500,900))
    

    plt_x = plot(t_fine_sim, Xsim', labels=[L"c_a" L"c_d" L"c_t"], title="States", xlabel="", ylabel="CO2 liq. conc. [kmol/m3]", ylims = (-1,5))
    plt_F = plot(t_vec, vcat(Usim[1,1], Usim[1,:]), seriestype = :steppre, title="Optimized Control Input", xlabel="Time (hr)", ylabel="Flow Rate", ylims=(0.2,0.45), label=L"F")
    plot!(plt_F, t_vec, vcat(lb[1], lb), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
    plot!(plt_F, t_vec, vcat(ub[1], ub), seriestype=:steppre, linestyle=:dash, color=:red, label = "")
    
    plt_Q = plot(t_vec, vcat(Usim[2,1], Usim[2,:]), seriestype = :steppre, label=L"Q", title="", xlabel="", ylabel="Reboiler duty [kW]", ylims = (20,35))
    plt_cgina = plot(t_vec, vcat(Dsim[1,1], Dsim[1,:]), seriestype = :steppre, label=L"c_{in,a}^g", title="Disturbances", xlabel="", ylabel="Flue gas conc. [mol/m3]", ylims=(4,5))
    plt_Fgina = plot(t_vec, vcat(Dsim[2,1], Dsim[2,:]), seriestype = :steppre, label=L"F_{in,a}^g", title="", xlabel="", ylabel="Liq. flow rate [m3/h]", ylims=(140,220))

    scale = 700.0
    H2W_ratio = 0.65

    plt_out = plot(plt_cap_eff, plt_F, plt_cgina, plt_x, plt_Q, plt_Fgina, 
        layout = (2,3), link=:x, 
        size = (scale, scale * H2W_ratio)
    )

    return plt_out


end


