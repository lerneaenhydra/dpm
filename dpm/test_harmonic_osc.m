clc
clear variables
format short eng

%Simple example of optimal control of an undamped harmonic oscillator (e.g.
%for a pendulum with angle and angular velocity states, limited torque
%input, bring pendulum to zero angle/zero speed as quickly as possible.

%DP solver set-up
dp_inp = dpm();

%Problem configuration
T_f = 3*pi;
%Contains the sample period
dp_inp.prb.T_s = pi/8;
%Contains the lower bound for all state variables
dp_inp.prb.X_l = [-1; -1];
%Contains the upper bound for all state variables
dp_inp.prb.X_h = [1; 1];

%The number of state variables
dp_inp.prb.N_x = 2;
%The number of control variables
dp_inp.prb.N_u = 1;
%Contains the number of time-steps to simulate
dp_inp.prb.N_t = round(T_f / dp_inp.prb.T_s)+1;

%Contains the lower bound for the terminal condition for states
dp_inp.prb.XT_l = [-inf; -inf];
%Contains the upper bound for the terminal condition for states
dp_inp.prb.XT_h = [inf; inf];

%Contains the upper and lower bounds for the initial condition for states
dp_inp.prb.X0_l = [0.5; 0.5];
%Contains the upper bound for the initial condition for states
dp_inp.prb.X0_h = [inf; inf];
%The number of grid points to generate for each state variable at each
%sample
dp_inp.prb.N_x_grid = [51; 51];
%The number of grid points to generate for each control variable at each
%sample
dp_inp.prb.N_u_grid = 101;


%Contains the lower bound for all inputs
dp_inp.prb.U_l = -1;
%Contains the upper bound for all inputs
dp_inp.prb.U_h = 1;

%The amount to scale the grid extent in each iteration, centered about the
%previous optimal path
dp_inp.sol.mu_grid_dec.x = 0.75;
dp_inp.sol.mu_grid_dec.u = 0.85;
dp_inp.sol.mu_grid_inc.x = 1.051;
dp_inp.sol.mu_grid_inc.u = 1.051;
%Termination threshold for maximum number of iterations
dp_inp.sol.iter_max = 15;
%Set to true to allow re-gridding the state variables after each iteration
dp_inp.sol.regrid_x = true;
%Set to true to allow re-gridding the control variables after each
%iteration
dp_inp.sol.regrid_u = true;
%Set true to enable debug mode (break execution on error/'unexpected' state
dp_inp.sol.debug = false;
%System configuration
dp_inp.sol.fun = @test_model_harmonic_osc;
dp_inp.sol.plotfun = @plot_iter;
dp_inp.sol.unique_thrs = 1.01;

%Interpolation mode to use. Set to a string, whose valid values depend on
%the chosen value of N_x as follows;
%	1D; All methods supported by interp1
%	>=2D; All methods supported by the griddedinterp class
dp_inp.sol.interpmode = 'linear';
%Extrapolation mode to use. Set similarly as the interpmode field.
dp_inp.sol.extrapmode = 'nearest';

dp_inp.sol.pen_norm = 'squaredeuclidean';
dp_inp.sol.pen_thrs = 2.^2;
dp_inp.sol.pen_fun_s = @(x) 1;
dp_inp.sol.pen_fun_a = @(x) 1;

%Model set-up for two-dimensional problem
A = [0, 1; -1, 0];	%Set up A matrix for an undamped harmonic oscillator
mod_consts.matexp = expm(A*dp_inp.prb.T_s);
mod_consts.t_f = dp_inp.prb.T_s;

h_plot = figure(1);
[res, grd, t, c, map] = dpm(dp_inp, mod_consts, h_plot);

figure(2);
subplot(1,2,1);
plot(t, reshape([res{1}.x],2,[]));
grid on;
title('Optimal state trajectory for initial and final grid');
xlabel('t');
ylabel('x');
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t, reshape([res{end}.x],2,[]), '--');
init_str = arrayfun(@(x) sprintf('initial state %d', x), 1:2, 'un', false);
opt_str = arrayfun(@(x) sprintf('optimal state %d', x), 1:2, 'un', false);
legend([init_str, opt_str]);
hold off;

subplot(1,2,2);
stairs(t, [res{1}.u]);
title('Optimal control signal for initial and final grid');
grid on;
xlabel('t');
ylabel('u');
hold on;
stairs(t, [res{end}.u], '--');
legend('u_{init}', 'u_{final}');
hold off;


if dp_inp.sol.iter_max > 1
	figure(3);
	dummy = c/min(c)-1;
	dummy(dummy == 0) = min(dummy(dummy ~= 0));
	semilogy(dummy);
	grid on;
	title('Cost relative to lowest determined cost');
	ylabel('Relative cost');
	xlabel('Iteration');
end

h_mapplot1 = figure(4);
y_str = 'State 2 value';
x_str = 'State 1 value';
plot_map(map{1}, dp_inp, 'h', h_mapplot1, 'n_arrows', 15, 'res', res{1}, ...
	'colorbar', true, ...
	'c_interp', 5, ...
	'bound_tightness', 0.5, ...
	'title', 'Optimal control map, first iteration', ...
	'xlabel', x_str, ...
	'ylabel', y_str, ...
	'clabel', 'Cumulative cost to end');

h_mapplot2 = figure(5);
plot_map(map{end}, dp_inp, 'h', h_mapplot2, 'n_arrows', 15, 'res', res{end}, ...
	'colorbar', true, ...
	'c_interp', 5, ...
	'bound_tightness', 0.5, ...
	'title', 'Optimal control map, last iteration', ...
	'xlabel', x_str, ...
	'ylabel', y_str, ...
	'clabel', 'Cumulative cost to end');


h_mapplot3 = figure(6);
plot_map(map{1}, dp_inp, 'h', h_mapplot3, 'n_arrows', 0, 'res', res{1}, ...
	'colorbar', true, ...
	'c_interp', 5, ...
	'bound_tightness', 0.5, ...
	'title', 'n_{oss}/n_{tot}, first iteration', ...
	'xlabel', x_str, ...
	'ylabel', y_str, ...
	'opt_uniqueness', true, ...
	'clabel', 'n_{oss}/n_{tot}');