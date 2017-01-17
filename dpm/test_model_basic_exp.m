function [x_nn1, x_nn2, c] = test_model_basic_exp(x_n1, x_n2, u_n1, u_n2, t, t_s)
% [x_nn1, x_nn2, c] = TEST_MODEL_BASIC_EXP(x_n1, x_n2, u_n1, u_n2, t, inp)
% Calculates the state resulting from applying inputs u_n1 and u_n2 to the
% system while at state x_n1 and x_n2 in a format that is compatible with
% GPU calculation (i.e. only using numerical types, each state/control
% variable as a seperate input, and only scalar optinal parameters.
%
%	Inputs: 
%	x_n1	State 1 at time t_n. Column array where each element
%			corresponds to the current value of state 1.
%	x_n2	State 2 at time t_n. Otherwise identical to x_n1. u_n1	Control
%			1 at time t_n. Column array where each elemenent corresponds to
%			control variable 1.
%	u_n2	Control 2 at time t_n. Otherwise identical to u_n1.
%	t		The current time
%
%	Additional parameters. NOTE: These must be scalars ordered in
%	alphabetical order!
%	t_s		The model sample rate
%
%
%	Returns:
%	x_nn1	State 1 value at time t_{n+1} given the current states and
%			controls. Should be of size [size(x_n1, 1), 1].
%	x_nn2	Equivalent to x_nn1 but for state 2.
%			
%	c		The cost required to move from 

%In this example, let the model be defined as [y1' = u1; y2' = u2] and let
%the cost function be defined as abs(y1) + abs(y2). (I.e. two completely
%decoupled systems).

x_nn1 = x_n1 + u_n1 * t_s;
x_nn2 = x_n2 + u_n2 * t_s;

c = (abs(x_n1) + abs(x_n2)) * t_s;

%Waste some CPU time to highlight the difference between CPU and GPU
%execution. Note that the use of a for loop biases the results to the CPUs
%benefit as branches are relatively more expensive on the GPU compared to
%the CPU.

foo = x_nn1;
for i=1:1e3
	foo = exp(foo);
end

end

