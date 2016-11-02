function [res_basic, grd_basic, t_basic, c_basic, map_basic, dp_inp_basic] = test_basic(varargin)
%Bare-bones example for a one-dimensional problem.

p = inputParser;
addParameter(p, 'pen_thrs', []);
parse(p, varargin{:});


if(isempty(varargin))
	clc
	clear variables
	format short eng	
	pen_thrs = 0;
else
	pen_thrs = p.Results.pen_thrs;
end


%DP solver set-up
dp_inp_basic = dpm();

%Problem configuration
T_f = 4;
%Contains the sample period
dp_inp_basic.prb.T_s = 0.5;
%Contains the lower bound for all state variables
dp_inp_basic.prb.X_l = -1;
%Contains the upper bound for all state variables
dp_inp_basic.prb.X_h = 2;

%The number of state variables
dp_inp_basic.prb.N_x = 1;
%The number of control variables
dp_inp_basic.prb.N_u = 1;
%Contains the number of time-steps to simulate
dp_inp_basic.prb.N_t = round(T_f / dp_inp_basic.prb.T_s)+1;


%Contains the lower bound for the terminal condition for states
dp_inp_basic.prb.XT_l = 1;
%Contains the upper bound for the terminal condition for states
dp_inp_basic.prb.XT_h = inf;

%Contains the upper and lower bounds for the initial condition for states
dp_inp_basic.prb.X0_l = 1;
%Contains the upper bound for the initial condition for states
dp_inp_basic.prb.X0_h = inf;

%Contains the lower bound for all inputs
dp_inp_basic.prb.U_l = -1;
%Contains the upper bound for all inputs
dp_inp_basic.prb.U_h = 1;

%The number of grid points to generate for each state variable at each
%sample
dp_inp_basic.prb.N_x_grid = 25;
%The number of grid points to generate for each control variable at each
%sample
dp_inp_basic.prb.N_u_grid = 15;

%The amount to scale the grid extent in each iteration, centered about the
%previous optimal path
dp_inp_basic.sol.mu_grid_dec = 0.75;
dp_inp_basic.sol.mu_grid_inc = 1.051;
%Termination threshold for maximum number of iterations
dp_inp_basic.sol.iter_max = 1;
%Set to true to allow re-gridding the state variables after each iteration
dp_inp_basic.sol.regrid_x = true;
%Set to true to allow re-gridding the control variables after each
%iteration
dp_inp_basic.sol.regrid_u = true;
%Set true to enable debug mode (break execution on error/'unexpected' state
dp_inp_basic.sol.debug = false;
%System configuration
dp_inp_basic.sol.fun = @test_model_1d;
dp_inp_basic.sol.plotfun = @plot_iter;

%Interpolation mode to use. Set to a string, whose valid values depend on
%the chosen value of N_x as follows;
%	1D; All methods supported by interp1
%	>=2D; All methods supported by the griddedinterp class
dp_inp_basic.sol.interpmode = 'linear';

dp_inp_basic.sol.pen_norm = 'squaredeuclidean';
dp_inp_basic.sol.pen_thrs = pen_thrs;
dp_inp_basic.sol.pen_fun_s = @(x) 1;
dp_inp_basic.sol.pen_fun_a = @(x) 1;

mod_consts.T_s = dp_inp_basic.prb.T_s;

h_plot = figure(1);
[res_basic, grd_basic, t_basic, c_basic, map_basic] = dpm(dp_inp_basic, mod_consts, h_plot);

figure(2);
subplot(1,2,1);
plot(t_basic, reshape([res_basic{1}.x],1,[]));
grid on;
title('Optimal state trajectory for initial and final grid');
xlabel('t');
ylabel('x');
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t_basic, reshape([res_basic{end}.x],1,[]), '--');
init_str = arrayfun(@(x) sprintf('initial state %d', x), 1, 'un', false);
opt_str = arrayfun(@(x) sprintf('optimal state %d', x), 1, 'un', false);
legend([init_str, opt_str]);
hold off;

subplot(1,2,2);
stairs(t_basic, [res_basic{1}.u]);
title('Optimal control signal for initial and final grid');
grid on;
xlabel('t');
ylabel('u');
hold on;
stairs(t_basic, [res_basic{end}.u], '--');
legend('u_{init}', 'u_{final}');
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
	'xlabel', 'State value', ...
	'ylabel', 'Time', ...
	'clabel', 'Cumulative cost to end');

if dp_inp_basic.sol.iter_max >= 5
	h_mapplot_final = figure(5);
	plot_map(map_basic{5}, dp_inp_basic, 'h', h_mapplot_final, 'n_arrows', 25, 'res', res_basic{5}, ...
		'colorbar', true, ...
		'c_interp', 3, ...
		'bound_tightness', 0.75, ...
		'title', 'Optimal control map for fifth iteration', ...
		'xlabel', 'State value', ...
		'ylabel', 'Time', ...
		'clabel', 'Cumulative cost to end');
end