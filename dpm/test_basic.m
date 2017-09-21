function [res_basic, grd_basic, t_basic, c_basic, map_basic, dp_inp_basic] = test_basic(varargin)
%Bare-bones example for the dpm solver using a two-state/two-control
%problem optionally calculated using the local parallel pool or on the GPU.
%See test_model_1d_exp for the dynamic model definition. See test_loose for
%an example of how to first solve a lower-dimensional problem and use its
%results to reduce the search range for applicable problems.
%To fairly compare the computational time for GPU calculations use the
%profiler and compare the time needed to execute line
% '[scat_x_nn_g, scat_c_g] = arrayfun(inp.sol.fun_exp, x_g, u_g, t_g, opts{:});'
%with
% '[scat_x_nn, scat_c] = inp.sol.fun(x_n, u_n, N.t(n_t), mod_consts);'
%or
% '[scat_x_nn_k, scat_c_k] = inp.sol.fun(x_n(idx_call,:), u_n(idx_call,:), N.t(n_t), mod_consts);'
%in dpm.m

p = inputParser;
addParameter(p, 'pen_thrs', 2.1, @(x) isscalar(x) && isfinite(x) && x >= 0);	%Override the default penalty threshold
addParameter(p, 'gpu_calc', false, @(x) isscalar(x) && islogical(x));			%Enable model GPU offloading. See the dpm help section for details on when this is beneficial.
addParameter(p, 'cpu_parallel', false, @(x) isscalar(x) && islogical(x));		%Enable model and interpolation parallelization on CPU. See the dpm help section for details on when this is beneficial.
addParameter(p, 'time_inv', false, @(x) isscalar(x) && islogical(x));			%Enable the model time invariance assumption. See the dpm help section for details on when this is beneficial.
parse(p, varargin{:});

%DP solver set-up
dp_inp_basic = dpm();

%Problem configuration
T_f = 4;
%Contains the sample period
dp_inp_basic.prb.T_s = 0.5;
%Contains the lower bound for all state variables
dp_inp_basic.prb.X_l = [-1; -1];
%Contains the upper bound for all state variables
dp_inp_basic.prb.X_h = [2; 2];

%The number of state variables
dp_inp_basic.prb.N_x = 2;
%The number of control variables
dp_inp_basic.prb.N_u = 2;
%Contains the number of time-steps to simulate
dp_inp_basic.prb.N_t = round(T_f / dp_inp_basic.prb.T_s)+1;


%Contains the lower bound for the terminal condition for states
dp_inp_basic.prb.XT_l = [1; 1];
%Contains the upper bound for the terminal condition for states
dp_inp_basic.prb.XT_h = [inf; inf];

%Contains the upper and lower bounds for the initial condition for states
dp_inp_basic.prb.X0_l = [1; 1];
%Contains the upper bound for the initial condition for states
dp_inp_basic.prb.X0_h = [inf; inf];

%Contains the lower bound for all inputs
dp_inp_basic.prb.U_l = [-1; -2];
%Contains the upper bound for all inputs
dp_inp_basic.prb.U_h = [1; 2];

%The number of grid points to generate for each state variable at each
%sample
dp_inp_basic.prb.N_x_grid = [25; 25];
%The number of grid points to generate for each control variable at each
%sample
dp_inp_basic.prb.N_u_grid = [25; 25];

%The amount to scale the grid extent in each iteration, centered about the
%previous optimal 5path
dp_inp_basic.sol.mu_grid_dec = 0.75;
dp_inp_basic.sol.mu_grid_inc = 1.051;
%Termination threshold for maximum number of iterations
dp_inp_basic.sol.iter_max = 30;
%Set to true to allow re-gridding the state variables after each iteration
dp_inp_basic.sol.regrid_x = true;
%Set to true to allow re-gridding the control variables after each
%iteration
dp_inp_basic.sol.regrid_u = true;
%Set true to enable debug mode (break execution on error/'unexpected' state
dp_inp_basic.sol.debug = false;
%System configuration
dp_inp_basic.sol.fun = @test_model_basic;
%Arbitrarily limit the number of state/control combinations to a value that
%ensures that multiple passes are required. This allows for a more
%effective illustration of the effectivness of the cpu_parallel
%functionality.
dp_inp_basic.sol.fun_maxcombs = min(1e5, round(prod(dp_inp_basic.prb.N_x_grid)*prod(dp_inp_basic.prb.N_u_grid)/50));
dp_inp_basic.sol.plotfun = @plot_iter;
dp_inp_basic.sol.time_inv = p.Results.time_inv;
dp_inp_basic.sol.cpu_parallel = p.Results.cpu_parallel;
%GPU configuration
dp_inp_basic.sol.gpu_enable = p.Results.gpu_calc;
dp_inp_basic.sol.gpu_enter = @single;	%Use single-precision floating-point variables for good performance on standard "gaming" GPUs
dp_inp_basic.sol.gpu_exit = @double;	%Convert data back to the default double-precision floats
dp_inp_basic.sol.fun_exp = @test_model_basic_exp;


%If using the model-time-invariance assumption the grid must be identical
%for all samples. This means that IDP only ever makes sense for
%applications where the states/controls are relatively constant compared to
%the initial state. To illustrate this, artificially inflate the initial
%search space by a factor of 100 (which ensures that the final result gives
%a trajectory that is fairly constant in this range).
if(p.Results.time_inv)
	dp_inp_basic.prb.X_l = dp_inp_basic.prb.X_l * 100;
	dp_inp_basic.prb.X_h = dp_inp_basic.prb.X_h * 100;
	dp_inp_basic.prb.U_l = dp_inp_basic.prb.U_l * 100;
	dp_inp_basic.prb.U_h = dp_inp_basic.prb.U_h * 100;
end

%Interpolation mode to use. Set to a string, whose valid values depend on
%the chosen value of N_x as follows;
%	1D; All methods supported by interp1
%	>=2D; All methods supported by the griddedinterp class
dp_inp_basic.sol.interpmode = 'linear';
%Extrapolation mode to use. Allowable values depend on N_x in the same way
%as for the interpolation mode
dp_inp_basic.sol.extrapmode = 'nearest';

dp_inp_basic.sol.pen_norm = 'squaredeuclidean';
dp_inp_basic.sol.pen_thrs = p.Results.pen_thrs;
dp_inp_basic.sol.pen_fun_s = @(x) 1;
dp_inp_basic.sol.pen_fun_a = @(x) 1;

mod_consts.T_s = dp_inp_basic.prb.T_s;

h_plot = figure(1);
[res_basic, grd_basic, t_basic, c_basic, map_basic] = dpm(dp_inp_basic, mod_consts, h_plot);

figure(2);
subplot(1,2,1);
plot(t_basic, reshape([res_basic{1}.x],dp_inp_basic.prb.N_x,[]));
grid on;
title('Optimal state trajectory for initial and final grid');
xlabel('t');
ylabel('x');
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t_basic, reshape([res_basic{end}.x],dp_inp_basic.prb.N_x,[]), '--');
init_str = arrayfun(@(x) sprintf('initial state %d', x), 1:dp_inp_basic.prb.N_x, 'un', false);
opt_str = arrayfun(@(x) sprintf('optimal state %d', x), 1:dp_inp_basic.prb.N_x, 'un', false);
legend([init_str, opt_str]);
hold off;

subplot(1,2,2);
plot(t_basic, reshape([res_basic{1}.u],dp_inp_basic.prb.N_u,[]));
title('Optimal control signal for initial and final grid');
grid on;
xlabel('t');
ylabel('u');
hold on;
plot(t_basic, reshape([res_basic{end}.u],dp_inp_basic.prb.N_u,[]), '--');
init_str = arrayfun(@(x) sprintf('initial control %d', x), 1:dp_inp_basic.prb.N_u, 'un', false);
opt_str = arrayfun(@(x) sprintf('optimal control %d', x), 1:dp_inp_basic.prb.N_u, 'un', false);
legend([init_str, opt_str]);
hold off;


if dp_inp_basic.sol.iter_max > 1
	figure(3);
	dummy = c_basic/min(c_basic)-1;
	dummy(dummy == 0) = min(dummy(dummy ~= 0));
	semilogy(dummy);
	grid on;
	title('Cost relative to analytic solution cost');
	ylabel('Relative cost');
	xlabel('Iteration');
end

h_mapplot_1 = figure(4);
plot_map(map_basic{1}, dp_inp_basic, 'h', h_mapplot_1, 'n_arrows', 25, 'res', res_basic{1}, ...
	'colorbar', true, ...
	'c_interp', 5, ...
	'bound_tightness', 0.75, ...
	'title', 'Optimal control map for first iteration', ...
	'xlabel', 'State 1', ...
	'ylabel', 'State 2', ...
	'clabel', 'Cumulative cost to end');

if dp_inp_basic.sol.iter_max >= 5
	h_mapplot_final = figure(5);
	plot_map(map_basic{5}, dp_inp_basic, 'h', h_mapplot_final, 'n_arrows', 25, 'res', res_basic{5}, ...
		'colorbar', true, ...
		'c_interp', 3, ...
		'bound_tightness', 0.75, ...
		'title', 'Optimal control map for fifth iteration', ...
		'xlabel', 'State 1', ...
		'ylabel', 'State 2', ...
		'clabel', 'Cumulative cost to end');
end