function success = cbc_picard_iteration_x(exp)
% CBC_PICARD_ITERATION_X  Perform a fixed-point (or Picard) iteration on
%   the higher harmonics of x.
%
%   SUCCESS = CBC_PICARD_ITERATION_X(EXP) returns true if the experiment
%   converged to a fixed point and false otherwise.
%
%   See also CBC_CREATE_EXPERIMENT.

% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015

% Make sure we are at steady-state
success = cbc_wait_for_convergence(exp);
if ~success
    warning('cbc_picard_iteration_x:convergence', '%s: Failed to converge to steady-state initially', datestr(now, 13));
    return;
end

% Get the current target
target = exp.rtc.par.x_target_coeffs;

% The harmonics of interest (include the DC level to remove any
% mis-alignment issues)
harmonics = [exp.fourier.idx_DC, exp.fourier.idx_higher];

% Do the Picard iteration
i = 0;
x = exp.rtc.par.x_coeffs_ave;
while (i==0) || ...
        (any(abs(target(harmonics) - x(harmonics)) > exp.opt.x_coeffs_tol) && (i < exp.opt.max_picard_iter))
    % Target is too far from the measured value (or on first iteration) so
    % update the target with the actual measured values
    target(harmonics) = x(harmonics);
    exp.rtc.par.x_target_coeffs = target;
    % Make sure we are at steady-state
    success = cbc_wait_for_convergence(exp);
    if ~success
        warning('picard_iteration_x:convergence', '%s: Failed to converge to steady-state', datestr(now, 13));
        return;
    end
    x = exp.rtc.par.x_coeffs_ave;
    i = i + 1;
end

% Did we converge in time?
if any(abs(target(harmonics) - x(harmonics)) > exp.opt.x_coeffs_tol)
    success = false;
else
    success = true;
end

end
