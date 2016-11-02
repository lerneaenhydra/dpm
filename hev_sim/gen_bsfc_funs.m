%Generate anonymous functions that calculate the BSFC and the
%instantaneuous fuel consumption for a given ICE.
%These functions are then used in the system model function as well as any
%other locations that need to calculate the BSFC and/or fuel consumption.
%Note that this method is used rather than an .m file with a function as
%the overhead for calling a function is significant.

ice_model = 'quadratic';	%Set to the model to use

keystr_bsfc = '%#!KEYSTR_BSFC_94ed5811c494d0e7186a6d05ab5e584e';
keystr_fuelinst = '%#!KEYSTR_FUELINST_fa8cce87b9df91aa087e5da136558638';

if(strcmp(ice_model, 'quadratic'))
	syms w_min k_w tau_min k_tau bsfc_min w tau
	bsfc = bsfc_min + sqrt((k_w * abs(w - w_min)).^2 + (k_tau * abs(tau - tau_min)).^2);
	
	bsfc = subs(bsfc, {w_min, k_w, tau_min, k_tau, bsfc_min}, {120, 2e-10, 50, 5e-10, gpkwh2kgpj(200)});
	
	fuelinst = simplify(bsfc * abs(w * tau));
	
	h_bsfc = matlabFunction(bsfc);
	h_fuelinst = matlabFunction(fuelinst);
else
	error('Invalid model type specified!')
end

str_bsfc = func2str(h_bsfc);
str_fuelinst = func2str(h_fuelinst);

str_bsfc = ['bsfc_f = ', str_bsfc, ';'];
str_fuelinst = ['fuelinst_f = ', str_fuelinst, ';'];

util_matched_replace('./', '.m', keystr_bsfc, str_bsfc);
util_matched_replace('./', '.m', keystr_fuelinst, str_fuelinst);