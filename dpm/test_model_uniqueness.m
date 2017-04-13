function [x_nn, c] = test_model_uniqueness(x_n, u_n, t, inp)
% [x_nn, c] = TEST_MODEL_BASIC(x_n, u_n, t, inp) Calculates the state
% resulting from applying a inputs u_n to the system while at states x_n.
% Calls the expanded system model for the actual computation.
%	x_n		State at time t_n. Array where columns correspond to state
%			variables. Each row corresponds to a given state configuration.
%	u_n		Control input at time t_n. Array where columns correspond to
%			control variables. Each row corresponds to a given control 
%			configuration.
%	t		The current time
%	inp		Optional input data.
%	x_nn	State values at time t_{n+1} given the state x_n at time t_n
%			and the control u_n. Apply calculations on a row-by-row basis.
%			(IE. each row in x_n and u_n should uniquely determine each row
%			in x_nn). Should be of size size(x_n).
%	c		The cost required to move from a given row in x_n to the same
%			row in x_nn. Set to inf if state and control results in an
%			invalid/undefined system transition.

%State variable is a simple integrator
x_nn = x_n + u_n * inp.prb.T_s;

%Apply a small constant cost so that we can compare all costs to each other
%(as we otherwise wouldn't be able to compare relative costs for any costs
%equal to zero).
c = abs(u_n) * inp.prb.T_s + realmin;

%Apply large cost if x_n is outside some forbidden zone 
out_bounds_l = inp.bounds_f(x_n, t);

out_bounds = zeros(size(out_bounds_l));
out_bounds(out_bounds_l) = inf;

c = c + out_bounds;

end

