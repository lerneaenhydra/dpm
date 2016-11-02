%GEN_BEST_FUELINST_OP Based on the current ICE model generate a function
%file that returns the minimum fuel consumption for a given demanded power.
%Also returns the motor's operating point (angular velocity and torque).

debug = false;

keystr = '%#!KEYSTR_OPT_FUELINST_0bb98bdd5464689e9ba6583ccc7ac141';

%Load constants
consts = param_setup('1D');

p_max = consts.shaft_maxw * consts.ice_maxtau;

p = linspace(0, p_max, consts.ice_fuelinst_opt_n);
w_best = zeros(size(p));
tau_best = zeros(size(p));

for n = 1:length(p)
	%Manually handle the case for zero power
	if(p(n) == 0)
		w_best(n) = 0;
		tau_best(n) = 0;
		continue;
	end
	
	%Find the minimum fuel consumption operating point for all the power demands
	mdot = @(w) calc_fuelinst(w, p(n)./w);
	%Set constraints so that we don't exceed the torque limit
	lb = p(n)./consts.ice_maxtau;
	ub = consts.shaft_maxw;
	x0 = mean([lb ub]);
	%Solve the constrained minimization. Use an inline function that scales
	%the fuel consumption to be 1 at the initial guess.
	opts = optimset('fmincon');
	opts.Display = 'none';
	w_best(n) = fmincon(@(x) mdot(x)*(1/mdot(x0)), x0, [], [], [], [], lb, ub, [], opts);
	tau_best(n) = p(n) / w_best(n);
	if debug
		clf;
		fprintf('Searching for best result with a demanded power of %f W\n', p(n));
		w_plot = linspace(lb, ub, 1e3);
		mdot_plot = mdot(w_plot);
		plot(w_plot, mdot_plot);
		hold on;
		best_mdot = mdot(w_best(n));
		scatter(w_best(n), best_mdot);
		title('Tested shaft velocity range and fuel consumption');
		xlabel('Angular velocity [rad/s]');
		ylabel('Fuel consumption [kg/s]');
		grid on;
		legend('Fuel consumption', 'Optimization result');
		keyboard;
	end
end

mdot_best = calc_fuelinst(w_best, tau_best);

p_str = mat2str(p);
w_str = mat2str(w_best);
mdot_str = mat2str(mdot_best);

replace_str = { ['p = ', p_str, ';'], ...
	['w = ', w_str, ';'], ...
	['mdot = ', mdot_str, ';'], ...
	'best_mdot = @(pwr) interp1(p, mdot, pwr, ''linear'', nan);', ...
	'best_w = @(pwr) interp1(p, w, pwr, ''linear'', nan);',...
};

util_matched_replace('./', '.m', keystr, replace_str);