function draw_sysplots(res_dpm, grd, c, t, iter, const, h)
%DRAW_SYSPLOTS Draws all system plots

figure(h);

t = t.';
res.st = cat(1, res_dpm.x);
res.ctrl = cat(1, res_dpm.u);
res.w = res.st(:,1);
res.soc = res.st(:,2);
res.tau_ice = res.ctrl(:,1);
res.tau_gen = res.ctrl(:,2);
res.tau_net = calc_tau_net(res.tau_ice, res.tau_gen, res.w, const);

subplot(3,2,1);
%Plot all torques
stairs(t, [res.tau_net, res.tau_ice, res.tau_gen]);
title('ICE/generator torques');
xlabel('time [s]');
ylabel('torque [Nm]');
grid on;
legend({'$\tau_{net}$','$\tau_{ICE}$','$\tau_{gen}$'},'Interpreter','latex');


subplot(3,2,2);
%Plot angular velocity
plot(t, res.w);
title('Generator/ICE angular velocity');
xlabel('time [s]');
ylabel('angular velocity [rad/s]');
grid on;

subplot(3,2,3);
%Plot BSFC chart and the system trajectory in it
hold off;
fcontour(@(x, y) kgpj2gpkwh(calc_bsfc(x, y)), [0, const.shaft_maxw, 0, const.ice_maxtau]);
colorbar;
title('ICE BSFC trajectory map [g/kWh]');
xlabel('Angular speed [rad/s]');
ylabel('Torque [Nm]');
grid on;
hold on;
%Place markers for the motor's operating points
scatter(res.w, res.tau_ice, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0.05);
hold off;

subplot(3,2,4);
%Plot the BSFC over time and the instantaneous power
yyaxis left;
plot(t, calc_fuelinst(res.w, res.tau_ice));
title('ICE BSFC and fuel consumption rate');
xlabel('Time [s]');
ylabel('Instantaneous fuel consumption [kg/s]');
grid on;

yyaxis right;
plot(t, kgpj2gpkwh(calc_bsfc(res.w, res.tau_ice)));
ylabel('BSFC [g/kWh]');


subplot(3,2,5);
%Plot all instantaneuous powers
res.pgen = -res.w .* res.tau_gen - res.tau_gen.^2 * const.gen_quadloss;
tot_p = const.tot_p_f(t);
tot_p = tot_p(:);
res.pbatt = tot_p - res.pgen;
plot(t, [res.pgen, tot_p, res.pbatt]);
title('Electrical subsystem instantaneuous power');
xlabel('time [s]');
ylabel('power [W]');
grid on;
legend({'$P_{gen}$','$P_{net}$','$P_{batt}$'}, 'Interpreter', 'latex');

subplot(3,2,6);
%Plot all battery data
yyaxis left;
plot(t, res.soc);
title('Battery state of charge and power');
xlabel('time [s]');
ylabel('soc [-]');
grid on;

yyaxis right;
stairs(t, res.pbatt);
ylabel('battery power [W]');
legend('SOC', 'net');

end

