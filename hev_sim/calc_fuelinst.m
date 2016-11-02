function fuel_inst = calc_fuelinst(w, tau)
%CALC_FUELINST Calculate the instantaneuous fuel consumption [kg/s]
%
%Function generated by gen_bsfc_funs.m
%DO NOT MODIFY THIS FUNCTION! Create function using gen_bsfc_funs.m!
%
%#!KEYSTR_FUELINST_fa8cce87b9df91aa087e5da136558638
fuelinst_f = @(tau,w)abs(tau.*w).*(sqrt(abs(tau-5.0e1).^2.*3.741444191567112e32+abs(w-1.2e2).^2.*5.986310706507379e31).*2.584939414228211e-26+5.555555555555555e-8);

fuel_inst = fuelinst_f(tau, w);

end

