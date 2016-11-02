function [x_nn, c] = test_model_1d(x_n, u_n, t, inp)
% [x_nn, c] = TEST_MODEL_1D(x_n, u_n, t, inp) Calculates the state resulting
% from applying an input u_n to the system while at state x_n.
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

%In this example, let the model be defined as y' = u and let the cost
%function be defined as abs(y).

x_nn = x_n + u_n * inp.T_s;

c = sum( abs(x_n) * inp.T_s, 2);


end

