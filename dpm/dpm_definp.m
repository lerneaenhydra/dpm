function [def_inp, grid_subset] = dpm_definp()
%[def_inp, grid_subset] = DPM_DEFINP() Returns an empty structure of the
%input used by the dpm function.
%	This function is internally called by the dpm script as needed.
% Copyright (c) 2016, Jonathan Lock
% All rights reserved.
%
% This file is part of DPM.
%
% DPM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% DPM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DPM.  If not, see <http://www.gnu.org/licenses/>.

%Contains the sample period
def_inp.prb.T_s = [];
%Contains the lower bound for all state variables
def_inp.prb.X_l = [];
%Contains the upper bound for all state variables
def_inp.prb.X_h = [];

%The number of state variables
def_inp.prb.N_x = [];
%The number of control variables
def_inp.prb.N_u = [];
%Contains the number of time-steps to simulate
def_inp.prb.N_t = [];


%Contains the lower bound for the terminal condition for states
def_inp.prb.XT_l = [];
%Contains the upper bound for the terminal condition for states
def_inp.prb.XT_h = [];

%Contains the upper and lower bounds for the initial condition for states
def_inp.prb.X0_l = [];
%Contains the upper bound for the initial condition for states
def_inp.prb.X0_h = [];

%Contains the lower bound for all inputs
def_inp.prb.U_l = [];
%Contains the upper bound for all inputs
def_inp.prb.U_h = [];

%The number of grid points to generate for each state variable at each
%sample
def_inp.prb.N_x_grid = [];
%The number of grid points to generate for each control variable at each
%sample
def_inp.prb.N_u_grid = [];

%Optionally manually force a number of state combinations to be included in
%the state grid. Should be a matrix with def_inp.prb.N_x rows and an
%arbitrary number of columns.
def_inp.prb.X_grid_manpts = [];
%Optionally manually force a number of control combinations to be included in
%the control grid. Should be a matrix with def_inp.prb.N_u rows and an
%arbitrary number of columns.
def_inp.prb.U_grid_manpts = [];

%Structure format to use for each element in def_inp.prb.grid_seed
%The type of variable (state or control) whose range is to be reduced. Set
%to 'x' for a state variable, or 'u' for a control variable.
grid_subset.vartype = [];
%Indexer for the specific variable to reduce the range for. E.g. set to 2
%to reduce the range for the second state/control variable.
grid_subset.varidx = [];
%The sample points (in time) that will be used to generate the non-uniform
%grid.
grid_subset.t = [];
%Center-points for grid generation for each time in grid_subset.t
grid_subset.center = [];
%Grid extent for grid generation for each time in grid_subset.t
grid_subset.range = [];
%Interpolation method used for grid generation between samples in
%grid_subset.t. May be any value that is accepted by the interp1 function.
grid_subset.interpmode = [];

%Optional non-uniform grid setup. Set to a cell array of grid_subset
%structures to force the dpm script to reduce the search space for any
%number of state/control variables. The order of the elements in the cell
%array has no effect.
def_inp.prb.grid_seed = cell(0);


%The amount to scale the grid extent in each iteration, centered about the
%previous optimal path depending on whether or not the previous iteration
%generated a valid solution or not. May be one of the following types;
%	scalar		- Scale all controls and states by the same amount
%	struct		- A struct with fields;
%					.x	A column vector with N_x scalars, where the n'th
%					element contains the amount to scale the n'th state
%					variable by
%					.u	A column vector with N_u scalars, where the n'th
%					element contains the amount to scale th n'th control by
def_inp.sol.mu_grid_dec = [];
def_inp.sol.mu_grid_inc = [];
%Termination threshold for maximum number of iterations
def_inp.sol.iter_max = [];
%Set to true to allow re-gridding the state variables after each iteration
def_inp.sol.regrid_x = [];
%Set to true to allow re-gridding the control variables after each
%iteration
def_inp.sol.regrid_u = [];

%Set true to enable debug mode (break execution on error/'unexpected'
%state)
def_inp.sol.debug = [];
%Set to any of;
%	'none'	Don't display any status messages
%	'final' Only keep the most recent iteration's status and cost
%	'iter'	Display each iteration's status and cost
def_inp.sol.display = [];

%System configuration
%The maximum number of state/control combinations to test per call to the
%system model. Reasonable values are typically on the order of 1e2 -- 1e9
def_inp.sol.fun_maxcombs = [];
%The system dynamics model. Should be a function of type
%[x_new, cost] = fun(x, u, t, opts) where;
%	x is an n by prb.N_x_grid array with current state values
%	u is an n by prb.N_u_grid array with control inputs
%	t is a scalar with the current time
%	opts is a struct with optional data
%	x_new is an n by prb.N_x_grid array with the state value at the next
%	sample after applying the control u to the current state x on a
%	row-by-row basis
%	cost is an n by 1 vector with the stage cost, i.e. the net cost of
%	applying u to x and arriving at x_new
%	Note: if using GPU calculation this function is typically only a
%	wrapper for the sol.fun_exp function. (See test_model_basic.m).
def_inp.sol.fun = [];
%A generic plot function, called once per IDP iteration
def_inp.sol.plotfun = [];

%GPU-related setup
%Set true to enable GPU model calculation
def_inp.sol.gpu_enable = [];
%Function to apply to input data before sending to GPU (e.g. @single)
def_inp.sol.gpu_enter = [];
%Function to apply to input data before sending to GPU (e.g. @double)
def_inp.sol.gpu_exit = [];
%Expanded system dynamics model. Should be a function of type
%[x_new, cost] = fun(x, u, t, a, b, c, ...) where the outputs and first
%three inputs are identical to sol.fun and a, b, c, ... are general
%numerical inputs (i.e. arrays). Typically these should be orderd
%alphabetically as sol.fun is simply a wrapper for this function.
def_inp.sol.fun_exp = [];

%Optional flag that, if present and set to true, indicates that the system
%model is time invariant (i.e. def_inp.sol.fun/def_inp.sol.fun_exp gives
%identical output for all tested input values t). If this is applicable and
%this flag is set the execution time will be increased very significantly.
%Note that this requires that the state/control grid is equal at all time
%instances, so this will give poor perfomance for applications where the
%state or control varies significantly over time, especially if used
%together with the range reducing method (i.e. non-empty
%def_inp.prb.grid_seed cells) or iterative dynamic programming method (i.e.
%def_inp.sol.iter_max set to any value other than 1).
def_inp.sol.time_inv = [];

%Optional vector that may be set to a sclar or vector of real-valued
%numbers typically > 0. If present will fill in the relative number of
%one-step-suboptimal controls that give a relative increase in cumulative
%cost less than each value listed in this scalar/vector. Here,
%one-step-suboptimal means applying a suboptimal control u to a given
%state x(n), generating a new state x(n+1), after which fully-optimal
%controls are taken. I.E. the returned value for each state (which, by
%nature, lies between 0 and 1) indicates how many of the tested controls
%give a small increase in the cumulative cost. This can give an idea as to
%if a problem is underdefined in certain regions/times. For example, if the
%returned values are always 0 this indicates that only the chosen control
%gives a cumulative cost that is within unique_thrs of the cumulative cost.
%However, if the returned value is instead, say, 0.5, this indicates that
%half of the tested control signals lie within unique_thrs of the
%cumulative cost and choosing a sub-optimal control would only result in
%small penalties to the cumulative cost.
%Typically this value is set to value(s) in the range 1e-6 to 1e-1. See
%'test_uniqueness.m' for an example of how to use this and view the results.
def_inp.sol.unique_thrs = [];

%Interpolation mode to use. Set to a string, whose valid values depend on
%the chosen value of N_x as follows;
%	1D; All methods supported by interp1
%	>=2D; All methods supported by the griddedinterp class
def_inp.sol.interpmode = [];

%Extrapolation mode to use, only used in forward-calculation phase. Set to
%a string, whose valid values depend on the chosen value of N_x as follows;
%	1D; All methods supported by interp1, typically 'inf' or 'extrap'
%	>=2D; All methods supported by the griddedinterp class, typically
%	'none' or 'nearest'
def_inp.sol.extrapmode = [];

%Norm to use for determining boundary for penalizing grid points near
%infeasible regions. Set to a string containing any of the norms supported
%by the pdist2 function.
def_inp.sol.pen_norm = [];
%Threshold to apply for penalty function; grid points further than this
%distance from the nearest infeasible point, as measured by the metric
%defined in def_inp.sol.pen_norm, will be completely unaffected by the
%penalization function.
def_inp.sol.pen_thrs = [];
%Penalization scaling function; the cumulative cost during the
%back-calculation phase for grid points closer than def_inp.sol.pen_thrs to
%any infeasible point will be scaled by this value, where the input to the
%function is the minimum distance to the nearest infeasible point. One
%example of a penalization function is
% @(dist) (def_inp.sol.pen_thrs - dist + 1)
%which will apply a linearly decreasing penalty that is equal to
%def_inp.sol.pen_thrs for feasible grid points that are neighbors with
%infeasible points.
def_inp.sol.pen_fun_s = [];
%Penalization additive function; functions similarly to pen_fun_s above,
%except adds a penalization term, rather than the multiplicative term 
%above.
def_inp.sol.pen_fun_a = [];

end

