function [x_nn, c] = test_model_1d_exp(x_n, u_n, t, t_s)
% [x_nn, c] = TEST_MODEL_1D_exp(x_n, u_n, t, inp) Calculates the state
% resulting from applying an input u_n to the system while at state x_n in
% a format that is compatible with GPU calculation (i.e. only using
% numerical types).
%
%	Inputs:
%	x_n		State at time t_n. Column array where each element corresponds
%			to the current state.
%	u_n		Control input at time t_n. Column array where each elemenent 
%			corresponds to the control variables.
%	t		The current time
%
%	Additional parameters. NOTE: These must be ordered in alphabetical
%	order!
%	t_s		The model sample rate
%
%
%	Returns:
%	x_nn	State values at time t_{n+1} given the state x_n at time t_n
%			and the control u_n. Apply calculations on a row-by-row basis.
%			(IE. each row in x_n and u_n should uniquely determine each row
%			in x_nn). Should be of size size(x_n).
%	c		The cost required to move from a given row in x_n to the same
%			row in x_nn. Set to inf if state and control results in an
%			invalid/undefined system transition.

%In this example, let the model be defined as y' = u and let the cost
%function be defined as abs(y).

x_nn = x_n + u_n * t_s;

c = abs(x_n) * t_s;

%Waste some CPU time to highlight the difference between CPU and GPU
%execution. Note that the use of a for loop biases the results to the CPUs
%benefit as branches are relatively more expensive on the GPU compared to
%the CPU.

foo = x_nn;
for i=1:1e3
	foo = exp(foo);
end

end

