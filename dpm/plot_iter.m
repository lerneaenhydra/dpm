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

	figure(h);
	
	clf;
	subplot(3,2,1);
	plotstuff(t, xo, state_finitespace);
	grid on;
	title('State progression and backwards-reachable state space');
	xlabel('Time');
	ylabel('State');

	subplot(3,2,2);
	plotstuff(t(1:end-1), uo, ctrl_range(1:end-1,:));
	grid on;
	title('Control');
	xlabel('Time');
	ylabel('Control');

	subplot(3,2,3);
	plot(t(1:end-1), c_samp, '-o');
	grid on;
	title('sample cost');
	xlabel('time');
	ylabel('cost');

	subplot(3,2,4);
	plot(t(1:end-1), cum, '-o');
	grid on;
	title('cumulative cost');
	xlabel('time');
	ylabel('cost');

	subplot(3,2,5:6);
	dummy = c(1:iter)/min(c(1:iter))-1;
	if(all(dummy == 0))
		dummy = ones(size(dummy));
	else
		dummy(dummy == 0) = min(dummy(dummy ~= 0));
	end
	semilogy(dummy, '-o');
	grid on;
	title('Cost relative to lowest determined cost');
	ylabel('Relative cost');
	xlabel('Iteration');


end

function plotstuff(t, var, bound)
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
		plot(t, var(:,k), '-o');
		th = plot(t(1:size(bound(:,k),1)), bound(:,k), 'LineStyle', '--');
		set(get(get(th,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
		th = plot(t(1:size(bound(:,k),1)), bound(:,k+N), 'LineStyle', '--');
		set(get(get(th,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	end
end