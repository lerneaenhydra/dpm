clc
format short eng

%In this script we will use a reduced search-space iterative dynamic
%programming method to (relatively speaking) quickly solve a loosely
%coupled system of two ODE's with wildly different time scales.
%We will see that this method converges to the (known) optimum despite
%using an initial approximation that neglects the dynamics of the ODE with
%faster time scale.
%
%For more details, see the IFAC conference paper "A Computationally Fast
%Iterative Dynamic Programming Method for Optimal Control of Loosely
%Coupled Dynamical Systems with Different Time Scales" Jonathan Lock, Tomas
%McKelvey, 2016.
%

%Define all project-wide constants
const.T_f = 3;				%Total simulation time
mod_const.alpha = 4;		%Time constant ratio between slow and fast integrators
const.T_s_1d = 0.5;			%Sample period for initial, 1-dimensional, problem. Must be a divisor of the total simulation time
const.T_s_2d = const.T_s_1d / mod_const.alpha;	%Sample period for full, 2-dimensional, problem. Must be a divisor of the total simulation time

if(rem(const.T_f,const.T_s_1d) ~= 0 || rem(const.T_f,const.T_s_2d) ~= 0 )
	error('Sample periods must be divisors of the total simulation time!');
end


%Set up the DP solver for the initial, 1-dimensional, case

inp = dpm();

inp.prb.T_s = const.T_s_1d;
inp.prb.X_l = -1;
inp.prb.X_h = 2;
inp.prb.N_x = 1;
inp.prb.N_u = 1;
inp.prb.N_t = const.T_f / const.T_s_1d + 1;
inp.prb.XT_l = 1;
inp.prb.XT_h = inf;
inp.prb.X0_l = 1;
inp.prb.X0_h = inf;
inp.prb.U_l = -1;
inp.prb.U_h = 1;
inp.prb.N_x_grid = 50;
inp.prb.N_u_grid = 50;

inp.sol.mu_grid_dec = 0.5;
inp.sol.mu_grid_inc = 1.051;
inp.sol.iter_max = 25;
inp.sol.regrid_x = true;
inp.sol.regrid_u = true;
inp.sol.debug = false;
inp.sol.fun = @model_loose_1d;
inp.sol.plotfun = @plot_iter;
inp.sol.interpmode = 'linear';
inp.sol.extrapmode = inf;
inp.sol.pen_norm = 'squaredeuclidean';
inp.sol.pen_thrs = sqrt(2)^2;
inp.sol.pen_fun_s = @(x) 2;
inp.sol.pen_fun_a = @(x) 1;

mod_const.T_s = const.T_s_1d;

h_iterplot = figure(1);

[res_1d, grd_1d, t_1d, c_1d, map_1d] = dpm(inp, mod_const, h_iterplot);

%Plot the deviation of the first and final iteration of the 1-dimensional
%case from the known (analytic) solution
an_sol = [linspace(inp.prb.X0_l, 0, abs(inp.prb.U_l)/inp.prb.T_s+1), ...
	zeros(1, inp.prb.N_t - 2*(abs(inp.prb.U_l)/inp.prb.T_s+1)), ...
	linspace(0, inp.prb.XT_l, abs(inp.prb.U_h)/inp.prb.T_s+1)];
figure(2);
hold off;
semilogy(t_1d, abs(an_sol - [res_1d{1}.x]));
grid on;
semilogy(t_1d, abs(an_sol - [res_1d{1}.x]));
hold on;
semilogy(t_1d, abs(an_sol - [res_1d{end}.x]), '--');
grid on;
str = sprintf('Deviation from analytic solution at first and %d''th iteration', inp.sol.iter_max);
title(str);
xlabel('Time');
ylabel('Deviation');
legend('First', 'Final', 'Location', 'Best');

if false
	h_mapplot = figure();
	plot_map(map_1d{1}, inp, 'bound_tightness', 0.5, 'n_arrows', inf, 'h', h_mapplot, 'res', res_1d{1});
end


%% Now, set up the DP solver for the full, 2-dimensional, case
%Initialize the search space based on the results from the previous,
%1-dimensional, problem

[inp, grid_subset] = dpm();

%
% To generate a valid solution we need to use penalized soft constraints
% for the fast state variable. This is equivalent to allowing a larger
% input magnitude u for the one-dimensional case, so one would expect the
% optimal solution to lie along these new, wider, bounds. However, this
% does not occur as we can apply a large, but finite, penalty to states
% outside the desired state range.
%
% We need to include this larger search space as we can only get a feasible
% solution if the 'fast' state variable exactly reaches +/- 1 (as for
% values whose magnitude is less than 1 we would quickly go outside the
% search space). This requirement means that we can only get feasible
% transitions if we can move from some state to exactly +/- 1, which
% requires us to be "lucky" in the sense that one of the gridded control
% values u "happens" to bring us to +/- 1. This is of course unlikely in
% the general sense which makes many of the states we would expect to be
% feasible infeasible.
%

inp.prb.T_s = const.T_s_2d;
inp.prb.X_l = [-1; -1.25];		%Note we allow values outside the nominal +/- 1 bound
inp.prb.X_h = [2; 1.25];
inp.prb.N_x = 2;
inp.prb.N_u = 1;
inp.prb.N_t = const.T_f / const.T_s_2d + 1;
inp.prb.XT_l = [1; -1];
inp.prb.XT_h = [inf; 0];
inp.prb.X0_l = [1; 0];
inp.prb.X0_h = [inf; 1];
inp.prb.U_l = -1;
inp.prb.U_h = 1;
inp.prb.N_x_grid = [50; 50];
inp.prb.N_u_grid = 25;

mod_const.T_s = const.T_s_2d;
mod_const.X_l_unpen = [-1; -1];	%Lower limit of state value range that is completely unpenalized
mod_const.X_h_unpen = [2; 1];	%Upper limit of state value range that is completely unpenalized
mod_const.pen = 1.1;				%Additive penalty to sample cost for violated soft constraints

%Feed the results from the previous, 1-dimensional problem, to the initial
%grid setup for this 2-dimensional problem. We want to reduce the search
%range for the first state variable based on res_1d{end}.x.
grid_subset.vartype = 'x';	%Set up the reduced grid range for a state variable
grid_subset.varidx = 1;		%Set the reduced grid range for the first variable
grid_subset.t = t_1d;		%The previous results were sampled at the previous sample rate
grid_subset.interpmode = 'linear';		%Apply linear interpolation to samples between the points in the previous solution
grid_subset.center = [res_1d{end}.x];	%The centerpoint for the variable search range is set to the previous solution
grid_subset.range = 1/mod_const.alpha * 2 * (inp.prb.X_h(1) - inp.prb.X_l(1)) * ones(size(t_1d));	%The variable range extent is set based on the time scale quotient alpha and the original search range
inp.prb.grid_seed{1} = grid_subset;		%Finally, set grid_seed to a 1x1 cell array with the configured variable and search range

inp.sol.mu_grid_dec = 0.8;
inp.sol.mu_grid_inc = 1.051;
inp.sol.iter_max = 50;
inp.sol.regrid_x = true;
inp.sol.regrid_u = true;
inp.sol.debug = false;
inp.sol.fun = @model_loose_2d;
inp.sol.plotfun = @plot_iter;
inp.sol.interpmode = 'linear';
inp.sol.extrapmode = 'none';

inp.sol.pen_norm = 'squaredeuclidean';
inp.sol.pen_thrs = sqrt(2)^2;
inp.sol.pen_fun_s = @(dist) 1;
%Apply an addiditive penalty to states close to the infeasibility bounds,
%otherwise we'll usually get an infeasible result during the
%fordward-calculation phase =(
inp.sol.pen_fun_a = @(dist) 1.1;

h_iterplot = figure(1);

[res_2d, grd_2d, t_2d, c_2d, map_2d] = dpm(inp, mod_const, h_iterplot);

h_mapplot = figure(2);
plot_map(map_2d{1}, inp, 'bound_tightness', 0.7, 'n_arrows', inf, 'h', h_mapplot, 'res', res_2d{1}, 'disp_grid', true);