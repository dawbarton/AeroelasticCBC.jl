norminf(x) = norm(x, Inf)
recordEQ(x, p) = (h=x.h,α=x.α,w=x.w)

function recordPO(x, p)
	xtt = BifurcationKit.getPeriodicOrbit(p.prob, x, p.p)
	period = BifurcationKit.getPeriod(p.prob, x, p.p)
	return (max = maximum(xtt[1,:]), min = minimum(xtt[1,:]), period = period)
end

function trace_from_equilibrium(p0)
    # Initial solution
    ax = Axis((h=1, h′=2, α=3, α′=4, w=5, w′=6))
    u0 = ComponentArray(zeros(6), ax)

    # Equilibrium problem
    prob = BifurcationProblem((u, p)->unsteady_model(ComponentArray(u, ax), p)[:], u0, p0, (@lens _.U);
        recordFromSolution = recordEQ)
    opts_br = ContinuationPar(pMin = 0.0, pMax = 25.0, ds = -0.1, dsmax = 0.1, nInversion = 8,
        detectBifurcation = 3, maxBisectionSteps = 25, nev = 3, plotEveryStep = 20,
        maxSteps = 1000, θ = 0.3)
    @set! opts_br.newtonOptions.verbose = false
    br_eq = continuation(prob, PALC(), opts_br; bothside = true, normC = norminf)

    # Periodic orbit problem
    optn_po = NewtonPar(verbose = false, tol = 1e-8,  maxIter = 25)
    opts_po_cont = ContinuationPar(dsmax = 0.1, ds= -0.01, dsmin = 1e-4, pMax = 25.0,
        pMin = 0.0, maxSteps = 200, newtonOptions = (@set optn_po.tol = 1e-8), nev = 3,
        tolStability = 1e-4, detectBifurcation = 3, plotEveryStep = 20, saveSolEveryStep=1,
        nInversion = 6)

    Mt = 50  # number of time sections
	br_po = continuation(
        br_eq, 2, opts_po_cont,
        PeriodicOrbitTrapProblem(M = Mt; updateSectionEveryStep = 1, jacobian = :Dense);
        ampfactor = 2., δp = 0.1,
        recordFromSolution = recordPO,
        finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
            # try to stop at the Hopf bifurcation to avoid continuing the trivial equilibrium
			return norm(z.u[begin:end-1]) > 1e-6
		end,
        normC = norminf, bothside = true)
    return (br_eq, br_po)
end
