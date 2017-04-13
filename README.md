# DPM(inp, mod_consts)
Simple dynamic programming solver implemented in Matlab for arbitrary continuous-valued problems represented with ODEs. Solves problems of the form:

min sum(c(x(k), u(k), k)) for k = [0, K]

subject to

* x(k+t) = f(x(k), u(k), k)
* b_in(x(k), u(k), k) <= 0
* x(k) is some m-dimensional vector
* u(k) is some l-dimensional vector

where x(k) and u(k) represent the sampled state and control variable
values at some time k.

This dynamic programming solver differs from the textbook DP solver in
the following ways;
	
* The search range for state and control variables is iteratively
 decreased following a mostly typical iterative dynamic programming
 (IDP) implementation. This solver differs from the typical IDP solver
 by allowing the grid size to increase slightly in the event that the
 feasible solution can be found after decreaseing the search grid
 extents.
	
* For loosely coupled multidimensional problems, an approximate
 solution generated by some n-1-dimensional problem can be used to
 reduce the search space for the n-dimensional problem, significantly
 reducing the time needed to arrive to a solution.
	
* An optional regularization term can by added and/or scale the cost
 function for tested states that are within some adjustable distance
 to infeasible state(s). This regularization term sigificanlty increases
 the ability for the DP solver to generate a feasible control
 trajectory for applications where the optimal control trajectory lies
 along a boundary of infeasibility (eg. bang-bang control).

Use the call `inp = dpm()` to return a structure containing the fields
required by the dpm solver.

Use the call `[inp, grid_subset] = dpm()` to return the previously
described structure as well as an empty structure of the form used to
configure limiting the search space for grid variable(s).

Use the call `[res, grid, t, c, map] = dpm(inp, mod_consts, h_iterplot)`
to solve the dynamic programming problem.

See the files `test_basic.m`, `test_loose.m`, `test_harmonic_osc.m`, and `test_uniqueness.m` for usage examples.
