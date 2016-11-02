function [x_nn, c] = model_loose_2d(x_n, u_n, t, inp)
%[x_nn, c] = MODEL_LOOSE_2D(x_n, u_n, t, inp)
%	Returns the next state x_nn and cost c when applying the control u to
%	the system state x_n at time t and with general inputs inp for the
%	full two-dimensional problem.
%
% The two-dimensional system is defined as
%
%	y2' = alpha * u;
%	y1' = y2;
%
%	with u constant during some interval [0, t_s]
%
%	where we want to solve u for
%
%	min(sum((|y1| + |y2/10|) * t_s, k = 0, k = N))
%
%	where t_s is the sampling period and subject to
%
%	y(0) = y0
%	y(N) = y_T
%	u_min <= u <= u_max
%
%	This set of ODE's can be analytically solved and the solution is;
%
%	y2(t) = y2(t_0) + alpha * t * u 
%	y1(t) = y1(t_0) + y2(t_0) * t + alpha * u * t^2 / 2
%
%	With this result it's easy to discretize our problem while getting the
%	same results as the continuous-time case.
%
%	To ensure the presence of a solution, we will allow y2 to go outside of
%	the bounds we will apply to the final problem. To ensure the result
%	converges to the bounds we wish to apply an additional penalty will be
%	applied to values that exceed the boundaries.
%
%

%Generate new state
x_nn = zeros(size(x_n));
x_nn(:,2) = x_n(:,2) + inp.alpha * inp.T_s * u_n;
x_nn(:,1) = x_n(:,1) + x_n(:,2) * inp.T_s + inp.alpha * u_n * inp.T_s.^2 / 2;

%Calculate base cost
c = sum((abs(x_n(:,1)) + abs(x_n(:,2))/10)* inp.T_s, 2);

%Penalize violation of soft constraints
idx_pen = any(x_nn < repmat(inp.X_l_unpen.', size(x_nn, 1),1), 2) | any(x_nn > repmat(inp.X_h_unpen.', size(x_nn, 1),1), 2);
c(idx_pen) = c(idx_pen) + inp.pen;


end