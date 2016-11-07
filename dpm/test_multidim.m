clc
clear variables
format short eng

test_dim = 2;		%Set to 1 or 2 to test a one or two (state) dimensional problem

test_dim = round(test_dim);

if(test_dim < 1 || test_dim > 2)
	error('invalid value for test_dim!');
end

%DP solver set-up
dp_inp = dpm();

%Problem configuration
if(test_dim == 1)
	T_f = 4;
	%Contains the sample period
	dp_inp.prb.T_s = 0.5;
	%Contains the lower bound for all state variables
	dp_inp.prb.X_l = -1 * ones(test_dim, 1);
	%Contains the upper bound for all state variables
	dp_inp.prb.X_h = 2 * ones(test_dim, 1);
elseif(test_dim == 2)
	T_f = 3*pi;
	%Contains the sample period
	dp_inp.prb.T_s = pi/8;
	%Contains the lower bound for all state variables
	dp_inp.prb.X_l = -1 * ones(test_dim, 1);
	%Contains the upper bound for all state variables
	dp_inp.prb.X_h = 1 * ones(test_dim, 1);
else
	error('invalid value for test_dim!');
end

%The number of state variables
dp_inp.prb.N_x = test_dim;
%The number of control variables
dp_inp.prb.N_u = 1;
%Contains the number of time-steps to simulate
dp_inp.prb.N_t = round(T_f / dp_inp.prb.T_s)+1;


%Contains the lower bound for the terminal condition for states
if(test_dim == 1)
	dp_inp.prb.XT_l = 1;
elseif(test_dim == 2)
	dp_inp.prb.XT_l = [-inf; -inf];
else
	error('Invalid value for test_dim!');
end
%Contains the upper bound for the terminal condition for states
dp_inp.prb.XT_h = inf * ones(test_dim, 1);

if(test_dim == 1)
	%Contains the upper and lower bounds for the initial condition for states
	dp_inp.prb.X0_l = 1;
	%The number of grid points to generate for each state variable at each
	%sample
	dp_inp.prb.N_x_grid = 21* ones(test_dim, 1);
	%The number of grid points to generate for each control variable at each
	%sample
	dp_inp.prb.N_u_grid = 15;
elseif(test_dim == 2)
	%Contains the upper and lower bounds for the initial condition for states
	dp_inp.prb.X0_l = [0.5; 0.5];
	%The number of grid points to generate for each state variable at each
	%sample
	dp_inp.prb.N_x_grid = 30 * ones(test_dim, 1);
	%The number of grid points to generate for each control variable at each
	%sample
	dp_inp.prb.N_u_grid = 3;
else
	error('Invalid value for test_dim!');
end
%Contains the upper bound for the initial condition for states
dp_inp.prb.X0_h = inf * ones(test_dim, 1);

%Contains the lower bound for all inputs
dp_inp.prb.U_l = -1;
%Contains the upper bound for all inputs
dp_inp.prb.U_h = 1;

%The amount to scale the grid extent in each iteration, centered about the
%previous optimal path
dp_inp.sol.mu_grid_dec = 0.75;
dp_inp.sol.mu_grid_inc = 1.051;
%Termination threshold for maximum number of iterations
dp_inp.sol.iter_max = 10;
%Set to true to allow re-gridding the state variables after each iteration
dp_inp.sol.regrid_x = true;
%Set to true to allow re-gridding the control variables after each
%iteration
dp_inp.sol.regrid_u = true;
%Set true to enable debug mode (break execution on error/'unexpected' state
dp_inp.sol.debug = false;
%System configuration
if(test_dim == 1)
	dp_inp.sol.fun = @test_model_1d;
elseif(test_dim == 2)
	dp_inp.sol.fun = @test_model_2d;
else
	error('Invalid value for test_dim!');
end
dp_inp.sol.plotfun = @plot_iter;

%Interpolation mode to use. Set to a string, whose valid values depend on
%the chosen value of N_x as follows;
%	1D; All methods supported by interp1
%	>=2D; All methods supported by the griddedinterp class
dp_inp.sol.interpmode = 'linear';
%Extrapolation mode to use. Set similarly as the intermode field.
if(test_dim == 1)
	dp_inp.sol.extrapmode = inf;
elseif(test_dim == 2)
	dp_inp.sol.extrapmode = 'none';
else
	error('Invalid value for test_dim!');
end

dp_inp.sol.pen_norm = 'squaredeuclidean';
dp_inp.sol.pen_thrs = 1.^2;
dp_inp.sol.pen_fun_s = @(x) 1;
dp_inp.sol.pen_fun_a = @(x) 1;

if(test_dim == 2)
	%Model set-up for two-dimensional problem
	A = [0, 1; -1, 0];	%Set up A matrix for an undamped harmonic oscillator
	mod_consts.matexp = expm(A*dp_inp.prb.T_s);
	mod_consts.t_f = dp_inp.prb.T_s;
elseif(test_dim == 1)
	%Model set-up for one-dimensional problem
	mod_consts.T_s = dp_inp.prb.T_s;
else
	error('Invalid value for test_dim!');
end

h_plot = figure(1);
[res, grd, t, c, map] = dpm(dp_inp, mod_consts, h_plot);

figure(2);
subplot(1,2,1);
plot(t, reshape([res{1}.x],test_dim,[]));
grid on;
title('Optimal state trajectory for initial and final grid');
xlabel('t');
ylabel('x');
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t, reshape([res{end}.x],test_dim,[]), '--');
init_str = arrayfun(@(x) sprintf('initial state %d', x), 1:test_dim, 'un', false);
opt_str = arrayfun(@(x) sprintf('optimal state %d', x), 1:test_dim, 'un', false);
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

h_mapplot = figure(4);
if(test_dim == 2)
	y_str = 'State 2 value';
	x_str = 'State 1 value';
elseif(test_dim == 1)
	y_str = 'State value';
	x_str = 'Time';
else
	error('Invalid value for test_dim!');
end
[~, h_c] = plot_map(map{10}, dp_inp, 'h', h_mapplot, 'n_arrows', 15, 'res', res{10}, ...
	'colorbar', true, ...
	'c_interp', 5, ...
	'bound_tightness', 0.5, ...
	'title', 'Optimal control map', ...
	'xlabel', x_str, ...
	'ylabel', y_str, ...
	'clabel', 'Cumulative cost to end');