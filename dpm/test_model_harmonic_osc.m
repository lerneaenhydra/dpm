function [x_nn, c] = test_model_harmonic_osc(x_n, u_n, t, inp)
% [x_nn, c] = TEST_MODEL_HARMONIC_OSC(x_n, u_n, t, inp) Calculates the
% state resulting from applying an input u_n to the system while at state
% x_n.
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

%In this example, let us define the model as an arbitrary set of time
%invariant ODE's of the form z' = Az + Bu where z contains N_x elements.
%The analytic solution to this is then given (see Beta, p 209) as
%z(t) = e^(At) * z_0 + e^(At) * int(e^-At * Bu, t)
%however, as we will assume that u is kept constant during each sample
%u(t) is constant, meaning that this expression reduces to
%z(t) = e^(At) * z_0 + Bu * t^2/2
%As e^(At) is constant this is precalculated to improve performance

x_nn = zeros(size(x_n));
c = zeros(size(x_n,1), 1);

for i = 1:size(x_n, 1)
	z_n = x_n(i,:).';
	g_n = [0; u_n(i,:)];
	z_nn = inp.matexp * z_n + g_n * inp.t_f^2/2;
	x_nn(i,:) = z_nn.';
	c(i) = norm(z_nn) * inp.t_f;
end

end

