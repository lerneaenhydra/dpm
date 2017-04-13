function [h_p, h_c] = plot_map(map, dp_inp, varargin)
%PLOT_MAP Plots the state transition map for 1- and 2-dimensional problems
%	map		Map structure for the solved problem.
%	dp_inp	Configuration structure used for DPM script.
%
%	The following parameters can additionally be supplied to alter/add
%	behavior beyond the bare-bones defaults using parameter/value pairs as
%	follows;
%	h		Will draw to the figure referenced by this handle
%	n		If the problem is 2-dimensional, will display the map for time
%			sample n and not add a slider bar for changing the sample.
%	n_arrows Will limit the number of arrows in each axis to the closest
%			evenly divisible value this value using nearest-neighbor
%			decimation. (IE. if the x and y axes have 50 points each and a
%			value of 25 is specified every other point will have a drawn
%			array, while a value of 20 will result in drawing an arrow for
%			every third point). Set to 0 to disable all arrows (except
%			those drawn with the 'res' argument). Leave unspecified or set
%			to inf to display all arrows.
%	res		If supplied should be the result given by the dpm function (IE.
%			[res, grd, t, c, map] = dpm(...)), will plot the forward
%			calculation state trajectory. For 1-dimensional cases this will
%			be a line from t_0 to t_f, while for 2-dimensional cases this
%			will be a short line segment for each time sample.
%	colorbar If true will add a color bar to the plot.
%	title	Will write a string to the axes title
%	xlabel	Will write a string to the axes x-label
%	ylabel	Will write a string to the axes y-label
%	clabel	Will write a string to the axes color bar
%	c_interp Positive scalar that will increase the detail level of data
%			used to draw the contour plot.
%	bound_tightness Scalar in range 0<= ... <= 1 that controls the
%			tightness of the contour plot interpolation, where smaller
%			values will generate larger regions and larger values will
%			generate smaller regions.
%	disp_grid Logical scalar, if supplied and true will display all grid
%			points.
%	h_p		Handle to the plot
%	h_c		If a color bar was added to the plot will be a handle to the
%			color bar.
%	opt_uniqueness Integer k, will override typical behavior and instead of
%			coloring the cumulative cost from each state to the terminal
%			state draw the relative number of suboptimal controls that,
%			if taken and followed by optimal controls, would not raise the total cost by more than 1 + dp_inp.sol.unique_thrs(k). Gives an
%			indication of the uniqueness of the problem
%			formulation, where a larger value for small thresholds
%			unique_thrs(k) indicates that variations around the generated
%			solution would typically give a similar final cost, while
%			smaller values indicate that the particular chosen optimal path
%			is significantly better than the other candidates.
%
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

	p = inputParser;
	addRequired(p, 'map');
	addRequired(p, 'dp_inp');
	addParameter(p, 'h', [], @(x) ishandle(x));
	addParameter(p, 'n', [], @(x) isnumeric(x) && isscalar(x) && round(x) == x);
	addParameter(p, 'n_arrows', inf, @(x) isnumeric(x) && isscalar(x) && round(x) == x && x >= 0);
	addParameter(p, 'res', [], @(x) isstruct(x));
	addParameter(p, 'colorbar', false, @(x) islogical(x) && isscalar(x));
	addParameter(p, 'title', '', @(x) ischar(x));
	addParameter(p, 'xlabel', '', @(x) ischar(x));
	addParameter(p, 'ylabel', '', @(x) ischar(x));
	addParameter(p, 'clabel', '', @(x) ischar(x));
	addParameter(p, 'c_interp', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
	addParameter(p, 'bound_tightness', 0.5, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
	addParameter(p, 'disp_grid', false, @(x) islogical(x) && isscalar(x));
	addParameter(p, 'opt_uniqueness', 0, @(x) x == round(x) && isscalar(x) && x > 0);
	
	parse(p, map, dp_inp, varargin{:});
	
	if(isempty(p.Results.h))
		h = figure();
	else
		h = p.Results.h;
		figure(h);
	end
	
	n = p.Results.n;
	n_arrows = p.Results.n_arrows;
	res = p.Results.res;
	en_colorbar = p.Results.colorbar;
	c_interp = p.Results.c_interp;
	bound_tightness = p.Results.bound_tightness;
	disp_grid = p.Results.disp_grid;
	opt_uniqueness = p.Results.opt_uniqueness;
	
	plotstr.title = p.Results.title;
	plotstr.xlabel = p.Results.xlabel;
	plotstr.ylabel = p.Results.ylabel;
	plotstr.clabel = p.Results.clabel;
	
	
	% Initialize plot and slider
	dim = dp_inp.prb.N_x;

	if(dim < 1 || dim > 2)
		error('Unsupported dimensionality!');
	end

	if(dim == 2)
		if(isempty(n))
			%Set up slider for the 2-dimensional case where there's no specific
			%sample requested
			
			%Clear current figure as the arrow plot function won't rescale the plot
			%axes
			clf;
			
			slider_max = size(map, 1);
			idx = 1;
			b = uicontrol('Parent', h, ...
			'Style', 'slider', ...
			'Units', 'Normalized', ...
			'Position', [0.05, 0.05, 0.02, 0.9], ...
			'value', idx, ...
			'SliderStep', [1/(slider_max-1), 10/(slider_max-1)], ...
			'min', 1, ...
			'max', slider_max, ...
			'Callback', @sliderCallback);
		else
			%If we've specified a certain sample to plot, do that
			idx = n;
		end
	elseif(dim == 1)
		%For the 1-dimensional case all samples will be plotted
		%simultaneously (ie. a y/t diagram) so the notion of an index
		%doesn't exist, however, set the idx variable to 1 to allow use of
		%the same code for both cases
		idx = 1;
	else
		error('unsupported dimensionality!');
	end

	%Set up slider callback to redraw plot on slider movement
	function sliderCallback(~, ~)
		figure(h);
		idx = round(get(b, 'value'));
		h_c = draw_map(map, dp_inp, idx, dim, n_arrows, res, en_colorbar, plotstr, c_interp, bound_tightness, disp_grid, opt_uniqueness);
	end

	%Force drawing plot on the first call to this function
	h_c = draw_map(map, dp_inp, idx, dim, n_arrows, res, en_colorbar, plotstr, c_interp, bound_tightness, disp_grid, opt_uniqueness);
	
	h_p = h;
end

function h_c = draw_map(map, dp_inp, idx, dim, n_arrows, res, en_colorbar, plotstr, c_interp, bound_tightness, disp_grid, opt_uniqueness)
	cla;
	hold off;
	%Generate the set of current/next state vectors to draw in the
	%arrow plot.
	if(n_arrows == inf || n_arrows == 0)
		arrow_decimation = ones(size(dp_inp.prb.N_x_grid.'));
	else
		arrow_decimation = ceil(dp_inp.prb.N_x_grid.' ./  n_arrows);
	end
	

	if(dim == 2)
		%Generate basic current/next state vectors and the associated
		%cumulative cost
		if(opt_uniqueness)
			contourbg_raw = map(idx).rel_unique_thrs(:,opt_uniqueness);
			cumcost = map(idx).cum(:);
			contourbg_raw(~isfinite(cumcost)) = inf;
		else
			contourbg_raw = map(idx).cum(:);
		end
		x = map(idx).x;
		xn = map(idx).xnn_scat;
		
		%Reshape the current state vectors to the format required by
		%the contourf function family
		x_c = reshape(x(:,1), dp_inp.prb.N_x_grid.');
		y_c = reshape(x(:,2), dp_inp.prb.N_x_grid.');
		c_c = reshape(contourbg_raw, dp_inp.prb.N_x_grid.');
		
		x_arrows = zeros([prod(ceil(dp_inp.prb.N_x_grid.'./arrow_decimation)), dim]);
		xn_arrows = zeros([prod(ceil(dp_inp.prb.N_x_grid.'./arrow_decimation)), dim]);
		
		for i=1:dim
			%Reshape the raw data to allow for easier decimation
			x_tmp = reshape(x(:,i), dp_inp.prb.N_x_grid.');
			xn_tmp = reshape(xn(:,i), dp_inp.prb.N_x_grid.');
			c_tmp = reshape(contourbg_raw, dp_inp.prb.N_x_grid.');
			
			%Decimate the raw data
			x_tmp = x_tmp(1:arrow_decimation(1):end, 1:arrow_decimation(2):end);
			xn_tmp = xn_tmp(1:arrow_decimation(1):end, 1:arrow_decimation(2):end);
			c_tmp = c_tmp(1:arrow_decimation(1):end, 1:arrow_decimation(2):end);
			
			%Reformat the decimated data and store in the final vector
			x_arrows(:,i) = x_tmp(:);
			xn_arrows(:,i) = xn_tmp(:);
			c_arrows = c_tmp(:);
		end
		
		%Finally, exclude the set of current/next state vectors with infinite cost
		x_arrows = x_arrows(isfinite(c_arrows), :);
		xn_arrows = xn_arrows(isfinite(c_arrows), :);
		
	elseif(dim == 1)
		%Extract raw data to plot
		if(opt_uniqueness)
			tmp = arrayfun(@(x) x.rel_unique_thrs(:,opt_uniqueness), map, 'un', false);
			contourbg_raw = cell2mat(tmp.');
			contourbg_raw(~isfinite([map.cum])) = inf;
		else
			contourbg_raw = [map.cum];
		end
		x = [map.x];
		xn = [map.xnn_scat];
		
		%For the 1-D case, the raw data is already in the format needed by
		%the contourf family, we just need to generate an array with time
		%data. As we want to display results up to and including the last
		%sample, extend the state and cumulative cost memory elements
		c_c = [contourbg_raw, zeros(size(contourbg_raw,1), 1)];
		y_c = [x, x(:,end)];
		t_c = linspace(0, dp_inp.prb.T_s * (dp_inp.prb.N_t-1), dp_inp.prb.N_t);
		x_c = repmat(t_c, size(y_c, 1), 1);
		
		%Generate data points for arrow diagram
		y_arrows = x(1:arrow_decimation:end, :);
		yn_arrows = xn(1:arrow_decimation:end, :);
		x_arrows = x_c(1:arrow_decimation:end, 1:end-1);
		xn_arrows = x_c(1:arrow_decimation:end, 2:end);
		
		c_arrows = contourbg_raw(1:arrow_decimation:end, :);
		
		x_arrows = [x_arrows(:), y_arrows(:)];
		xn_arrows = [xn_arrows(:), yn_arrows(:)];
		
		x_arrows = x_arrows(isfinite(c_arrows(:)), :);
		xn_arrows = xn_arrows(isfinite(c_arrows(:)), :);
		
		%For the contour plot we don't really have anything useful to
		%display during the last time sample, so remove these elements from
		%the data passed to the contourf function
		x_c = x_c(:,1:end-1);
		y_c = y_c(:,1:end-1);
		c_c = c_c(:,1:end-1);
	else
		error('unsupported dimensionality');
	end

	x_c_finite = x_c(isfinite(c_c));
	y_c_finite = y_c(isfinite(c_c));
	c_c_finite = c_c(isfinite(c_c));
	
	%Generate colormap to generate a mild gradient for cost and
	%indicate infeasible states and zero-cost states with a distinct color

	%Number of colors bands to draw for contour plot
	bands = 100;
	cmap = repmat(linspace(1, 0.5, bands).', 1, 3);
	
	%Set true if a contour map was drawn (which also sets the axis extents)
	drew_contour = false;

	%Draw contour map if there are enough finite costs to display, IE. we
	%can at least generate a polygon with nonzero area
	if(sum(sum(isfinite(c_c_finite))) >= 3)
		%The contourf function doesn't interpolate values in a particularly
		%pretty way, so expand the x_c, y_c, and c_c matrices so that the
		%effects are less prominent in the plot. Apply linear extrapolation as
		%we'll manually draw bounds around the feasible set later on, so we
		%don't want to highlight regions near the boundary yet.

		interp_c = scatteredInterpolant(x_c_finite, y_c_finite, c_c_finite, 'linear', 'nearest');

		x_c_interp = linspace(min(min(x_c)), max(max(x_c)), round(size(x_c, 2) * c_interp));
		x_c_interp = repmat(x_c_interp, round(size(x_c, 1) * c_interp), 1);

		y_c_interp = linspace(min(min(y_c)), max(max(y_c)), round(size(y_c, 1) * c_interp)).';
		y_c_interp = repmat(y_c_interp, 1, round(size(y_c, 2) * c_interp));

		warning('off', 'all');
		c_c_interp = interp_c(x_c_interp, y_c_interp);
		warning('on', 'all');
		%Check if the interpolation succeeded, this will fail if all points
		%in x_c_interp and y_c_interp are colinear. If we failed don't
		%attempt to draw a contour plot
		if(~isempty(c_c_interp))
			c_c_interp(~isfinite(c_c_interp)) = inf;
			contour_range  = linspace(min(min(c_c_interp(isfinite(c_c_interp)))), max(max(c_c_interp(isfinite(c_c_interp)))), bands);
			hold on;
			contourf(x_c_interp, y_c_interp, c_c_interp, contour_range, 'edgecolor', 'none');
			drew_contour = true;
			if(opt_uniqueness)
				colormap(flipud(parula));
			else
				colormap(cmap);
			end
			if(en_colorbar)
				h_c = colorbar;
				ylabel(h_c, plotstr.clabel);
			else
				h_c = [];
			end
		else
			h_c = [];
		end
	else
		h_c = [];
	end

	
	%Now manually generate the boundary for the unreachable/infeasible
	%space and superimpose this on the contour plot.
	k_b = boundary(x_c_finite, y_c_finite, bound_tightness);
	
	x_bound = x_c_finite(k_b);
	y_bound = y_c_finite(k_b);
	
	%As we want the region generated by x_bound and y_bound to be a hole we
	%need to create a rectangle with vertices at the current grid extents.
	%As the map at the current index contains all tested states we can
	%get the grid extents by simply looking at the minimum and maximum
	%values found in each column in x.
	if(dim == 1)
		x_lim = xlim;
		y_lim = ylim;
	elseif(dim == 2)
		x_lim = [min(x(:,1)), max(x(:,1))];
		y_lim = [min(x(:,2)), max(x(:,2))];
	else
		error('Unsupported dimensionality');
	end
	x_rect = [x_lim, fliplr(x_lim), x_lim(1)];
	y_rect = [repmat(y_lim(1), 1, 2), repmat(y_lim(2), 1, 2), y_lim(1)];
	
	[x_rect, y_rect] = poly2cw(x_rect, y_rect);
	[x_bound, y_bound] = poly2ccw(x_bound, y_bound);
	
	[f, v] = poly2fv({x_bound, x_rect}, {y_bound, y_rect});
	
	patch('Faces', f, 'Vertices', v, 'FaceColor', [0.6, 0.4, 0.4], 'EdgeColor', 'none');
	if(~drew_contour || dim == 2)
		axis([x_lim, y_lim]);
	end
	hold on;
	
	
	%Draw arrow-plot if there are any arrows to add
	if(~isempty(x_arrows) && n_arrows ~= 0)
		warning('off', 'all');
		arrow(x_arrows, xn_arrows, 'color', [0, 0.2, 0.5], 'length', 10);
		warning('on', 'all');
	end
	
	%If desired, plot the forward-calculation state trajectory
	if(~isempty(res))
		if(dim == 1)
			dummy = [res.x].';
			%For the 1D-case, plot forward trajectory with solid
			%(non-arrow) line as this is (arguably) more clear.
			plot(t_c, dummy, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2);
		elseif(dim == 2)
			x_fw = res(idx).x;
			xn_fw = res(idx+1).x;
			if(~isempty(x_fw) && all(all(isfinite(x_fw))))
				warning('off', 'all');
				arrow(x_fw, xn_fw, 'color', [0.1, 0.6, 0.2], 'Width', 2, 'length', 10);
				warning('on', 'all');
			end
		else
			error('unsupported dimensionality');
		end
	end
	
	%If desired, add draw points corresponding to the grid points
	if(disp_grid)
		scatter(x(:,1), x(:,2), 10, 'filled', 'k', 'MarkerEdgeColor', 'k')
	end
	
	hold off;
	
	if(~isempty(plotstr.title))
		title(plotstr.title);
	end
	if(~isempty(plotstr.xlabel))
		xlabel(plotstr.xlabel);
	end
	if(~isempty(plotstr.ylabel))
		ylabel(plotstr.ylabel);
	end
	
	grid on;
	axis tight;
end

