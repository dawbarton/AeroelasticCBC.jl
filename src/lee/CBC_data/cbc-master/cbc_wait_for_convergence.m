function success = cbc_wait_for_convergence(exp)
% CBC_WAIT_FOR_CONVERGENCE  Wait for the (normalise) variance of the
%   Fourier coefficients to settle to the desired level.
%
%   success = CBC_WAIT_FOR_CONVERGENCE(EXP) returns true if the Fourier
%   coefficients settled in the time allowed. 
%
%   See also CBC_CREATE_EXPERIMENT.

% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015

% Always wait to allow changes to propogate through the system and effect
% the variance calculations
pause(exp.opt.wait_time);

% Get the handle to the RTC
rtc = exp.rtc;

% Wait until the absolute variance or the normalised variance is
% sufficiently small (either being true is sufficient)
iter = 1;
x_coeffs_ave_norm = norm(rtc.par.x_coeffs_ave(rtc.fourier.idx_AC)); % Don't use the offset for normalisation
x_coeffs_var = rtc.par.x_coeffs_var;
while (iter < exp.opt.max_waits) ...
        && any(x_coeffs_var > exp.opt.x_coeffs_var_tol_abs) ...
        && any(x_coeffs_var/x_coeffs_ave_norm > exp.opt.x_coeffs_var_tol_rel)
    pause(exp.opt.wait_time);
    iter = iter + 1;
    x_coeffs_ave_norm = norm(rtc.par.x_coeffs_ave(rtc.fourier.idx_AC)); % Don't use the offset for normalisation
    x_coeffs_var = rtc.par.x_coeffs_var;
end

if any(x_coeffs_var > exp.opt.x_coeffs_var_tol_abs) ...
        && any(x_coeffs_var/x_coeffs_ave_norm > exp.opt.x_coeffs_var_tol_rel)
    success = false;
else
    success = true;
end

end
