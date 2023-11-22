export figXX_bifn_diagram

function figXX_bifn_diagram()
    p0 = unsteady_model_case2()
    # Extract the bifurcations from the bifurcation diagram
    (br_eq, br_po) = trace_from_equilibrium(p0)
    hopf = br_eq.specialpoint[2]
    (hopf.type == :hopf) || error("Not a Hopf bifurcation as expected - parameters changed?")
    sn = br_po.specialpoint[4]
    (sn.type == :bp) || error("Not a saddle-node bifurcation as expected - parameters changed?")
    U_hopf = hopf.param
    h_hopf = 0.0
    U_po = [br_po[i].param for i in 1:length(br_po)]
    h_po = [br_po[i].max for i in 1:length(br_po)]
    U_sn = sn.param
    h_sn = sn.printsol.max
    # Extract two example periodic orbits
    (br_eq, br_po) = trace_from_equilibrium(p0; save_par = 18.0)
    pt1 = br_po.specialpoint[2]
    (pt1.type == :var"save-1") || error("Not a periodic orbit as expected - parameters changed?")
    pt2 = br_po.specialpoint[4]
    (pt2.type == :var"save-1") || error("Not a periodic orbit as expected - parameters changed?")
    U_pt1 = pt1.param
    h_pt1 = pt1.printsol.max
    U_pt2 = pt2.param
    h_pt2 = pt2.printsol.max
    # Plot the bifurcation diagram

end
