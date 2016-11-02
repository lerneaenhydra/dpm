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

clear variables
close all
format short eng
%%
%Add the path to the dpm function
addpath([pwd filesep '..' filesep 'dpm' filesep]);


%Note; Use the gen_battmodel.m, gen_bsfc_funs.m, gen_best_fuelinst_op.m and
%gen_bsfc_map.m scripts to programmatically update the rest of the files in
%the workspace. Inlined anonymous functions are used in several locations
%rather than a function in an m-file as the performance penalty for calling
%an external m-function can be prohibitive otherwise.

%First, solve the 1-dimensional case. This is used to generate an
%approximate SOC trajectory that will then be used to generate the initial
%grid for the 2-dimensional problem where crankshaft dynamics are added.
consts_1d = param_setup('1D');

%Generate handles to used figures
h_iterplot_1d = figure(1);

[res_1d, grd_1d, t_1d, c_1d, map_1d] = dpm(consts_1d.dp_inp, consts_1d, h_iterplot_1d);
%%
save('1d_res.mat');
%%
consts_2d = param_setup('2D', res_1d{end});

h_iterplot_2d = figure(2);

[res_2d, grd_2d, t_2d, c_2d, map_2d] = dpm(consts_2d.dp_inp, consts_2d, h_iterplot_2d);

h_resplot_2d = figure(3);
draw_sysplots(res_2d{1}, grd_2d, c_2d, t_2d, 0, consts_2d, h_resplot_2d);

%%
%To reduce disk space usage, only save the first and last system map and
%grid
for i = 2:(length(map_2d)-1)
	map_2d{i} = [];
	grd_2d{i} = [];
end
save('2d_res.mat');