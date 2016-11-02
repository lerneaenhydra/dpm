%CALC_POWER_CYCLE Updates calc_power_cycle.m to contain the desired power
%cycle.
%Generated using
%http://www.virtual-car.org/wheels/wheels-road-load-calculation.html.

keystr = '%#!KEYSTR_POWER_CYCLE_cf618f1716decc2fa8b7c5e47d300b78';

%Read raw data and convert to meaningful data
rawdata = csvread('power_cycle_rawdata_US06.txt', 1, 0);

t_cyc = rawdata(:,1);
t_cyc = t_cyc - t_cyc(1);
p_neg = rawdata(:,9);
p_pos = rawdata(:,10);
p_cyc = p_pos + p_neg;

replace_str = {['t_cyc = ', mat2str(t_cyc), ';'], ...
	['p_cyc = ', mat2str(p_cyc), ';'], ...
	'p = interp1(t_cyc, p_cyc, t, ''linear'', nan);'};

util_matched_replace('./','.m',keystr,replace_str);

fprintf('Selected cycle has a test period of %f seconds. Be sure to update the total simulation time to match this!\n', t_cyc(end));