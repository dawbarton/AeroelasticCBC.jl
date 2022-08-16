function all_plots()
end

function plot_theme()
    return @pgf {
                 enlarge_x_limits = false,
                 font = "\\sffamily"
                 }
end

function figX_timeseries()
    p0 = unsteady_model_case2()
    p0.U = 20
    u0 = [0.1, 0, 0, 0, 0, 0]
    prob = ODEProblem(unsteady_model!, u0, (0.0, 100.0), p0)
    sol0 = solve(prob, Tsit5(); save_everystep = false, save_end = true)
    sol1 = solve(remake(prob; tspan = (0.0, 1.0), u0 = sol0.u[end]), Tsit5())
end
