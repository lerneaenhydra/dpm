%Simple test script to display the control uniqueness for a simple problem.
%Assume the problem to solve is the optimal control of a dynamic system of
%form x(n+1) = x(n) + u(n)*T_s. Require x(t) > a for t1 < t < t2 and x(t) <
%a for t3 < t < t4. Solely penalize abs(u). Finally display a plot
%indicating, for every state value x(n) and sample n, the number of
%one-step suboptimal controls that can be applied that would result in a
%total increase in cost less than some threshold. Here, one-step suboptimal
%implies taking one suboptimal control to x(n+1), and from there following
%the optimal control sequence to the terminal sample.

clc
clear variables
format short eng

dp_inp = dpm();

%Set up system for a sample rate of 0.1s, over a total time of 10s, with
%a search range of -1 <= x <= 1 and -1 <= u <= 1. Don't apply any
%initial/terminal constraints.

t_f = 10;

dp_inp.prb.T_s = 0.1;
dp_inp.prb.X_l = -1;
dp_inp.prb.X_h = 1;
dp_inp.prb.N_x = 1;
dp_inp.prb.N_u = 1;
dp_inp.prb.N_t = ceil(t_f / dp_inp.prb.T_s);
dp_inp.prb.XT_l = -1;
dp_inp.prb.XT_h = 1;
dp_inp.prb.X0_l = -1;
dp_inp.prb.X0_h = 1;
dp_inp.prb.U_l = -1;
dp_inp.prb.U_h = 1;
dp_inp.prb.N_x_grid = 101;
dp_inp.prb.N_u_grid = 1001;

dp_inp.sol.iter_max = 1;
dp_inp.sol.fun = @test_model_uniqueness;
dp_inp.sol.plotfun = @plot_iter;
dp_inp.sol.unique_thrs = 1e-3;
dp_inp.sol.interpmode = 'linear';
dp_inp.sol.extrapmode = nan;
dp_inp.sol.pen_norm = 'squaredeuclidean';
dp_inp.sol.pen_thrs = 0;			%Do not apply any regularization to states near infeasibility
dp_inp.sol.pen_fun_s = @(x) 1;
dp_inp.sol.pen_fun_a = @(x) 0;

h_iterplot = figure(1);
mod_consts = dp_inp;

a = 0.3;
t_bound = [2,3,7,8];
mod_consts.bounds_f = @(x, t) (x < a & t > t_bound(1) & t < t_bound(2)) | (x > -a & t > t_bound(3) & t < t_bound(4));

[res, grid, t, c, map] = dpm(dp_inp, mod_consts, h_iterplot);

%The first plot displays the minimum cumulative cost to go from any state
%value at any time to the final time, along with the optimal trajectory in
%green and infeasible regions in red.

%The second plot displays the same trajectory and infeasible regions,
%replacing the minimum cumulative cost with the relative number of
%suboptimal controls that, if taken, and then followed by taking from the
%resulting state's optimal controls, results in an increased cumulative
%cost of less than 'unique_thrs'.

%The second plot shows four distinct feasible regions;
%	1) for t >= 7 or where t>=3 and x <= -0.3 there only exists one
%	control that lies within 0.1% of the minimum cost, which is visible as
%	n_oss/n_tot == 0. This is fully reasonable, as there is no cost
%	associated with x and no terminal constraints are applied. Clearly, the
%	best control to apply is u=0, and anything else will be significantly
%	worse.
%	2) for all t <= 3 and x == 0.3, again, there only exists a single
%	control that is close-to optimal, as we require x >= 0.3 at 2 < t < 3,
%	so if we are already at x == 0.3 we do not want to decrease x (as we
%	would then need to increase it to at least 0.3 before t == 2) and we do
%	not want to increase it (as we require x <= -0.3 before t == 7, i.e. we
%	would need to decrease it again). Therefore, these is only one control
%	that is a good choice.
%	3) for all regions where n_oss/n_tot == 0.5 half of all the tested
%	controls have virtually the same cumulative cost. As the bounds on the
%	control signal are relatively large there is plenty of time to bring
%	the state variable to values that fulfils the constraints on x.
%	Furthermore, as the penalty is linear with u's magnitude it doesn't
%	matter how quickly the state is brought to values that fulfils the
%	constraints; so long as any control with the right sign is chosen the
%	cumulative cost will be the same (if ignoring the effects of
%	state/control quantization, which are small here). Because of this the
%	state trajectory for 3 < t < 7 is close-to underdetermined, and the
%	particular chosen trajectory is primarily controlled by the grid
%	quantization effects. (Which can be seen as the trajectory here changes
%	significantly when N_u and/or N_x is changed, but still remains
%	monotonically decreasing and always only just brings x to -0.3.)
%	4) for all t <= 3 and 0.2 <= x <= 0.4, and x ~= 0.3 there exists a
%	region with a decreasing number of controls that have close-to the same
%	cost. This is reasonable, as the number of tested controls that doesn't
%	overshoot x == 3 decreases linearly with abs(x-3).

plot_map(map{1}, dp_inp, 'res', res{1}, 'n_arrows', 0, 'bound_tightness', ...
	0.9, 'xlabel', 'Time [s]', 'ylabel', 'State value [-]', 'Title', ...
	'Cumulative cost and infeasible region map', 'colorbar', true, ...
	'clabel', 'Cumulative cost for optimal control');
for k = 1:length(dp_inp.sol.unique_thrs)
	titstr = sprintf('n_{oss} / t_{tot} s.t. C_{rel} < %g', dp_inp.sol.unique_thrs(k));
	plot_map(map{1}, dp_inp, 'opt_uniqueness', k, 'n_arrows', 0, ...
		'bound_tightness', 0.9, 'title', titstr, 'xlabel', 'Time [s]', ...
		'ylabel', 'State value [-]', 'colorbar', true, 'res', res{1});
end

