function tau_net = calc_tau_net(tau_ice, tau_gen, w, shaft_fric, shaft_fric_norm_w)
%CALC_TAU_NET Calculates the net torque for a given ICE and generator
%torque
tau_net = (tau_ice + tau_gen - shaft_fric * min(1, w/shaft_fric_norm_w) );
end

