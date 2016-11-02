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


%Structure format to use for each non-empty element in
%def_inp.prb.grid_seed
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

%Optional non-uniform grid setup. Set the n'th cell to a structure with the
%same fields as in grid_subset to force the dpm script to reduce the search
%space for the n'th state variable.
def_inp.prb.grid_seed = cell(0);


%The amount to scale the grid extent in each iteration, centered about the
%previous optimal path
def_inp.sol.mu_grid_dec = [];
def_inp.sol.mu_grid_inc = [];
%Termination threshold for maximum number of iterations
def_inp.sol.iter_max = [];
%Set to true to allow re-gridding the state variables after each iteration
def_inp.sol.regrid_x = [];
%Set to true to allow re-gridding the control variables after each
%iteration
def_inp.sol.regrid_u = [];
%Set true to enable debug mode (break execution on error/'unexpected' state
def_inp.sol.debug = [];
%System configuration
def_inp.sol.fun = [];
def_inp.sol.plotfun = [];

%Interpolation mode to use. Set to a string, whose valid values depend on
%the chosen value of N_x as follows;
%	1D; All methods supported by interp1
%	>=2D; All methods supported by the griddedinterp class
def_inp.sol.interpmode = [];

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

