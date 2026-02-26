# Applications

## What the Solver Handles

The solver finds elementary shortest paths on an undirected graph with a
single capacity resource and negative-cost edges.  The negative costs are
what makes the problem hard: they create negative-cost cycles, so LP
relaxations produce fractional subtours that must be cut off.  Without
negative costs the problem reduces to standard shortest-path and needs no
branch-and-cut.

The formulation uses vertex profits and vertex demands, but this is a
modelling convenience, not a structural limitation.


### Encoding Principles

**Vertex profits are negative-cost edges.**
Split a profitable node $v$ into $v_\text{in}$ and $v_\text{out}$ connected
by an edge with cost $-\pi_v$.  The profit is collected exactly when the
node is traversed.  The graph stays undirected and symmetric.

**Edge resource consumption becomes vertex demand.**
An edge-based budget constraint $\sum t_e x_e \le T$ can be moved onto
vertices by the same splitting: the internal edge $v_\text{in}$–$v_\text{out}$
carries the travel-time cost as a demand on $v$.  The capacity constraint
$\sum d_i y_i \le Q$ then enforces the budget.

These two observations mean that any problem with a symmetric graph, a
single resource, optional vertices, and potential negative-cost cycles can
be encoded into the solver's native formulation.

---

## Direct Applications

### Column Generation Pricing

The primary application.  In branch-and-price for vehicle routing, the
master LP produces dual values that become vertex profits in the pricing
subproblem.  Because duals can be large, negative-cost cycles appear
naturally, and the pricing solver must find **elementary** shortest paths.

Standard approaches use dynamic-programming labeling
(Feillet et al. 2004, Righini & Salani 2006).
This solver provides a branch-and-cut alternative following
Jepsen et al. (2014).

| Master problem | Pricing structure | Resource |
|---|---|---|
| CVRP (Capacitated VRP) | Symmetric, single capacity | Weight / demand |
| VRP with profits | Symmetric, single capacity | Weight / demand |

### Prize-Collecting TSP

Select a subset of customers to visit, collecting a prize at each, subject
to a penalty for skipped customers or a constraint on total prizes
collected.  The symmetric PCTSP (Balas 1989) maps directly: prizes are
profits, the penalty/constraint is the capacity resource.  Negative-cost
cycles arise whenever a cluster of nearby customers has prizes exceeding
the travel cost to visit them.

The undirected Selective TSP (Gendreau et al. 1998) is the same structure
with a different objective weighting.

### Orienteering Problem

A traveller with a distance or time budget visits locations to collect
scores.  The budget constraint $\sum t_e x_e \le T$ is an edge resource,
but via node splitting it encodes as vertex demands against a capacity $Q = T$.
Scores are vertex profits (negative-cost edges after splitting).

Negative-cost cycles appear when clustered high-value locations can be
visited for less travel cost than their combined score.  The branch-and-cut
approach with SEC and RCI cuts applies directly.

This covers a range of practical routing-with-profit problems
(Vansteenwegen et al. 2011, Archetti et al. 2014, Kobeaga et al. 2024):

- **Drone inspection** — battery budget, information value per site
- **Tourist trip planning** — time budget, attraction scores
- **Post-disaster assessment** — flight-time budget, damage-information value
- **Surveillance routing** — fuel budget, target reward

In all cases the distinguishing feature is that profitable clusters create
negative-cost cycles, so elementary-path enforcement via cutting planes is
needed.

---

## Extensions

Two extensions cover a large additional class of applications.
Neither requires fundamental changes to the solver architecture:
the MIP formulation, SEC separation, and branch-and-cut framework
carry over.  The main work is in adapting the capacity-sensitive
cuts (RCI, Multistar/GLM) and the labeling-based preprocessor.


### Directed Graphs

The undirected edge formulation (degree constraints) becomes a directed
flow formulation (flow conservation).  SEC separation via Gomory-Hu trees
generalises to directed min-cuts.  RCI and Multistar adapt to directed
demand flow.

This opens:

| Application | Why directed |
|---|---|
| Asymmetric PCTSP (Dell'Amico et al. 1995) | Sequence-dependent setup costs |
| VRPTW pricing | Directed arcs with time resources |
| Steel scheduling (Balas 1989) | Asymmetric changeover costs |
| General asymmetric RCSPP (Beasley & Christofides 1989) | Asymmetric networks |


### Multiple Resources

Adding a second capacity constraint $\sum d^2_i y_i \le Q_2$ is one
linear constraint in the MIP.  SEC cuts are resource-independent.
RCI gains a second rounding dimension (the Chvátal-Gomory argument
extends to each resource independently).  The labeling preprocessor
tracks an additional component per label.

This opens:

| Application | Resources |
|---|---|
| VRPTW pricing (capacity + time) | Weight, duration |
| Multi-compartment VRP pricing | Volume per compartment |
| EV routing with load | Weight, battery |
| Drone routing with payload | Payload weight, flight energy |

For VRPTW pricing in particular, the combination of directed graphs
and two resources (capacity and time) would cover the pricing
subproblem structure used by most modern exact VRP solvers
(Desaulniers et al. 2011, Pessoa et al. 2020).


### What We Do Not Target

**Time windows** require per-node interval constraints
($a_i \le \tau_i \le b_i$) that interact with the arrival time at each
node along the path.  This is inherently sequential and handled far more
naturally by dynamic-programming labeling (Feillet et al. 2004,
Righini & Salani 2006) than by cutting planes.  Branch-and-cut is the
wrong tool for time-window feasibility.

---

## References

- Archetti, C., Speranza, M. G., & Vigo, D. (2014). [Vehicle routing problems with profits](https://doi.org/10.1137/1.9781611973594.ch10). In *Vehicle Routing: Problems, Methods, and Applications* (2nd ed., pp. 273–298). SIAM.
- Balas, E. (1989). [The prize collecting traveling salesman problem](https://doi.org/10.1002/net.3230190602). *Networks*, 19(6), 621–636.
- Beasley, J. E. & Christofides, N. (1989). [An algorithm for the resource constrained shortest path problem](https://doi.org/10.1002/net.3230190402). *Networks*, 19(4), 379–394.
- Dell'Amico, M., Maffioli, F., & Värbrand, P. (1995). [On prize-collecting tours and the asymmetric travelling salesman problem](https://doi.org/10.1016/0969-6016(95)00010-5). *International Transactions in Operational Research*, 2(3), 297–308.
- Desaulniers, G., Desrosiers, J., & Spoorendonk, S. (2011). [Cutting planes for branch-and-price algorithms](https://doi.org/10.1002/net.20471). *Networks*, 58(4), 301–310.
- Feillet, D., Dejax, P., Gendreau, M., & Gueguen, C. (2004). [An exact algorithm for the elementary shortest path problem with resource constraints](https://doi.org/10.1002/net.20033). *Networks*, 44(3), 216–229.
- Gendreau, M., Laporte, G., & Semet, F. (1998). [A branch-and-cut algorithm for the undirected selective traveling salesman problem](https://doi.org/10.1002/(SICI)1097-0037(199812)32:4<263::AID-NET3>3.0.CO;2-Q). *Networks*, 32(4), 263–273.
- Jepsen, M., Petersen, B., Spoorendonk, S., & Pisinger, D. (2014). [A branch-and-cut algorithm for the capacitated profitable tour problem](https://doi.org/10.1016/S1572-5286(14)00036-X). *Discrete Optimization*, 14, 78–96.
- Kobeaga, G., Rojas-Delgado, J., Merino, M., & Lozano, J. A. (2024). [A revisited branch-and-cut algorithm for large-scale orienteering problems](https://doi.org/10.1016/j.ejor.2023.07.034). *European Journal of Operational Research*, 313(1), 44–68.
- Pessoa, A., Sadykov, R., Uchoa, E., & Vanderbeck, F. (2020). [A generic exact solver for vehicle routing and related problems](https://doi.org/10.1007/s10107-020-01523-z). *Mathematical Programming*, 183, 483–523.
- Righini, G. & Salani, M. (2006). [Symmetry helps: bounded bi-directional dynamic programming for the elementary shortest path problem with resource constraints](https://doi.org/10.1016/j.disopt.2006.05.007). *Discrete Optimization*, 3(3), 255–273.
- Vansteenwegen, P., Souffriau, W., & Van Oudheusden, D. (2011). [The orienteering problem: A survey](https://doi.org/10.1016/j.ejor.2010.03.045). *European Journal of Operational Research*, 209(1), 1–10.
