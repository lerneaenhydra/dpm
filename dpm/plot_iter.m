function plot_iter(res, grd, c, t, iter, const, h)
%PLOT_RES Generates plots based on the output from the DP solver
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

	xo = cat(1, res.x);
	uo = cat(1, res(1:end-1).u);
	c_samp = cat(1, res(1:end-1).c);
	cum = cat(1, res(1:end-1).cum);
	state_finitespace = cat(1, res.finitespace);
	ctrl_range = [cat(2,grd.u_h).', cat(2, grd.u_l).'];

	if(length(t > 50))
		plotstr = '-';
	else
		plotstr = '-o';
	end
	
	figure(h);
	
	clf;
	subplot(3,2,1);
	cla;
	plotstuff(t, xo, state_finitespace, plotstr);
	grid on;
	title('State progression and backwards-reachable state space');
	xlabel('Time');
	ylabel('State');
	axis auto;

	subplot(3,2,2);
	cla;
	plotstuff(t(1:end-1), uo, ctrl_range(1:end-1,:), plotstr);
	grid on;
	title('Optimal control and search space range');
	xlabel('Time');
	ylabel('Control');
	axis auto;

	subplot(3,2,3);
	cla;
	plot(t(1:end-1), c_samp, plotstr);
	grid on;
	title('sample cost');
	xlabel('time');
	ylabel('cost');
	axis auto;

	subplot(3,2,4);
	cla;
	plot(t(1:end-1), cum, plotstr);
	grid on;
	title('cumulative cost');
	xlabel('time');
	ylabel('cost');
	axis auto;

	subplot(3,2,5:6);
	cla;
	dummy = c(1:iter)/min(c(1:iter))-1;
	if(all(dummy == 0))
		dummy = ones(size(dummy));
	else
		dummy(dummy == 0) = min(dummy(dummy ~= 0));
	end
	if(length(c) > 50)
		plotstr = '-';
	else
		plotstr = '-o';
	end
	semilogy(dummy, plotstr);
	grid on;
	title('Cost relative to lowest determined cost');
	ylabel('Relative cost');
	xlabel('Iteration');
	axis auto;


end

function plotstuff(t, var, bound, plotstr)
	hold on;
	N = size(var, 2);
	colmap = parula(N + 1);
	for k = 1:N
		if N == 2
			if k == 1
				yyaxis left;
			else
				yyaxis right;
			end
		end
		plot(t, var(:,k), plotstr);
		th = plot(t(1:size(bound(:,k),1)), bound(:,k), 'LineStyle', '--');
		set(get(get(th,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
		th = plot(t(1:size(bound(:,k),1)), bound(:,k+N), 'LineStyle', '--');
		set(get(get(th,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	end
end