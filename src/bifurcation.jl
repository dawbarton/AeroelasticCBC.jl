norminf(x) = norm(x, Inf)
recordEQ(x, p) = (h = x[1], α = x[2], w = x[3])

function recordPO(x, p)
    xtt = BifurcationKit.get_periodic_orbit(p.prob, x, p.p)
    period = BifurcationKit.getperiod(p.prob, x, p.p)
    return (max = maximum(xtt[1, :]), min = minimum(xtt[1, :]), period = period,
            u0 = xtt[:, 1])
end

function trace_from_equilibrium(p0; save_par = nothing)
    # Initial solution
    u0 = zeros(6)

    # Equilibrium problem
    prob = BifurcationProblem(unsteady_model, u0, p0, (@lens _.U);
                              record_from_solution = recordEQ)
    opts_br = ContinuationPar(p_min = 0.0, p_max = 25.0, ds = -0.1, dsmax = 0.1,
                              n_inversion = 8,
                              detect_bifurcation = 3, max_bisection_steps = 25, nev = 3,
                              plot_every_step = 20,
                              max_steps = 1000, a = 0.3)
    @set! opts_br.newton_options.verbose = false
    br_eq = continuation(prob, PALC(), opts_br; bothside = true, normC = norminf)

    # Periodic orbit problem
    optn_po = NewtonPar(verbose = false, tol = 1e-8, max_iterations = 25)
    if save_par !== nothing
        detect_bifurcation = 0
        detect_event = 1  # only works with 1 (nearest neighbour); 2 (bisection) fails
        event = SaveAtEvent((save_par, 0.0))  # needs dummy parameter to prevent error
    else
        detect_bifurcation = 3
        detect_event = 0
        event = SaveAtEvent((0.0,))
    end
    opts_po_cont = ContinuationPar(dsmax = 0.1, ds = -0.01, dsmin = 1e-4, p_max = 25.0,
                                   p_min = 0.0, max_steps = 200,
                                   newton_options = (@set optn_po.tol = 1e-8), nev = 3,
                                   tol_stability = 1e-4, detect_bifurcation = detect_bifurcation,
                                   plot_every_step = 20, save_sol_every_step = 1,
                                   n_inversion = 6, detect_event = detect_event)

    Mt = 50  # number of time sections
    br_po = continuation(br_eq, 2, opts_po_cont,
                         PeriodicOrbitTrapProblem(M = Mt; update_section_every_step = 1,
                                                  jacobian = :Dense);
                         ampfactor = 2.0, δp = 0.1,
                         record_from_solution = recordPO,
                         finalise_solution = (z, tau, step, cont_result; prob = nothing, kwargs...) -> begin
                         # try to stop at the Hopf bifurcation to avoid continuing the trivial equilibrium
                         return norm(z.u[begin:(end - 1)]) > 1e-6 end,
                         normC = norminf, bothside = true, event = event)
    return (br_eq, br_po)
end
