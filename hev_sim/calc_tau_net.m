function tau_net = calc_tau_net(tau_ice, tau_gen, w, consts)
%CALC_TAU_NET Calculates the net torque for a given ICE and generator
%torque
tau_net = (tau_ice + tau_gen - consts.shaft_fric * min(1, w/consts.shaft_fric_norm_w) );
end

