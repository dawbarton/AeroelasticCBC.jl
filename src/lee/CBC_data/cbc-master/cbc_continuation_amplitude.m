function exp = cbc_continuation_amplitude(exp, varargin)
% CBC_CONTINUATION_AMPLITUDE  Do a continuation in forcing amplitude by
% continuously incrementing the amplitude of the control target.
%
% EXP = CBC_CONTINUATION_AMPLITUDE(RTC) performs the continuation using the
% default parameters and returns an experiment data structure containing
% the measured data.
%
% EXP = CBC_CONTINUATION_AMPLITUDE(EXP) performs the continuation using the
% parameters found in the experiment data structure and returns an
% experiment data structure containing the measured data.
%
% EXP = CBC_CONTINUATION_AMPLITUDE(..., Name, Value) performs the
% continuation as above with the named pairs of parameters and values
% overridden with the values supplied.
%
% Options
%
%     direction: allowed values are +1 or -1. Default +1.
%         Set the direction of the continuation.
%
%     out_range: allowed values are [a, b] where a < b. Default [0, Inf].
%         Set the range of allowable output amplitudes.
%
%     x_range: allowed values are [a, b] where a < b. Default [0, Inf].
%         Set the range of allowable x amplitudes.
%
%     adaptive: allowed values are true or false. Default true.
%         Use adaptive step sizes.
%
%     h0: allowed values are a > 0. Default 0.01.
%         Set the initial step size (if adaptive = false this is the fixed
%         step size).
%
%     h_max_scale: allowed values are a > 0. Default 1.5.
%         The maximum amount by which the step size can be scaled each
%         iteration (adaptive = true only).
%
%     step_range: allowed values are [a, b] where a < b. Default [0, 1].
%         The range bounding the norm of the step (adaptive = true only).
%
%     step_norm: allowed values are [a, b]. Default [1, 1].
%         The norm of the step is calculated as 
%           sqrt((step_norm(1)*delta_x)^2 + (step_norm(2)*delta_out)^2)
%         to allow for adjustment if the response and the input are on
%         different scales (adaptive = true only).
%
%     debug: allowed values are true or false. Default false.
%         Print the step size at each iteration (adaptive = true only).
%
%     plotting: allowed values are true or false. Default true.
%         Plot the amplitudes of the input and response as the continuation
%         runs.
%
%     hold: allowed values are true or false. Default false.
%         If hold is false, a new figure is created. If hold is true, the
%         current figure/axes is used.
%
%     colour: allowed values are valid Matlab colours. Default 'b'.
%         Set the plotting colour. Use in conjunction with hold.
%
% To disable step rejection warnings use 
%   warning('OFF', 'cbc_continuation_amplitude:rejected_step')
%
% See also CBC_CREATE_EXPERIMENT, WARNING.

% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015

% Parse the input
p = inputParser();
p.KeepUnmatched = true;
if ismethod(p, 'addParameter')
    % New versions of Matlab
    add_par = @p.addParameter;
else
    % Old versions of Matlab
    add_par = @p.addParamValue; %#ok<NVREPL>
end
p.addOptional('direction', 1, @(x)or(x == -1, x == +1));
add_par('out_range', [0, Inf], @(x)(x(2) > x(1)));
add_par('x_range', [0, Inf], @(x)(x(2) > x(1)));
add_par('h0', 0.01, @(x)(x > 0));
add_par('h_max_scale', 1.5, @(x)(x > 0));
add_par('step_range', [0, 1], @(x)(x(2) > x(1)));
add_par('step_norm', [1, 1]);
add_par('adaptive', true, @islogical);
add_par('plotting', true, @islogical);
add_par('hold', false, @islogical);
add_par('colour', 'b');
add_par('debug', false, @islogical);
add_par('max_fails', 5, @(x)(x > 0));
p.parse(varargin{:});

% Get the RTC handle
if ~isfield(exp, 'rtc')
    exp = cbc_create_experiment(exp, p.Unmatched); % Pass through unused options
else
    if ~isempty(fieldnames(p.Unmatched))
        unmatched = fieldnames(p.Unmatched)';
        unmatched(2, 1:end-1) = {' '};
        error('Unknown parameter(s): %s', [unmatched{:}]);
    end
end
rtc = exp.rtc;

% Check that the controller is active
if ~rtc.par.x_control
    error('Controller is not active');
end

% Set the name of the experiment
exp.name = sprintf([mfilename '_%04d%02d%02d_%02d%02d%02.0f'], exp.timestamp_created);

% Add the calling parameters to the experiment
exp.cont_par = p.Results;

% Get the indices of the fundamental harmonic
idx_sin = exp.fourier.idx_fund(1);
idx_cos = exp.fourier.idx_fund(2);

% Check that the target doesn't have any cosine component
if rtc.par.x_target_coeffs(idx_cos) ~= 0
    warning([mfilename ':target_coeffs'], 'Target has components in cosine as well as sine - did you forget to clear x_target_coeffs?');
end

% Check that we are converged to a steady-state
if ~cbc_wait_for_convergence(exp)
    error('%s: Failed to reach steady-state on first point!', datestr(now, 13));
end

% Get the current operating point
pt = cbc_get_data_point(exp);
pt.x_amp = norm(pt.x_coeffs_ave(exp.fourier.idx_AC));
pt.out_amp = norm(pt.out_coeffs_ave(exp.fourier.idx_AC));
exp.data = pt;

% Get the current target (sin coefficient only)
target = rtc.par.x_target_coeffs(idx_sin);
lasttarget = target;

% Initial step size
if exp.cont_par.direction >= 0
    h = exp.cont_par.h0;
else
    h = -exp.cont_par.h0;
end

% Ideal step
step_ideal = mean(exp.cont_par.step_range);

% Create a figure
if exp.cont_par.plotting
    if ~exp.cont_par.hold
        figure;
        ax = axes;
        xlabel('Forcing amplitude');
        ylabel('Response amplitude');
    else
        ax = gca;
    end
    hold(ax, 'on');
end

% Create a data directory
[s, msg] = mkdir('data');
if ~s
    error('Failed to create data directory; mkdir gave the error "%s"', msg);
end

% Iterate - 
% Check that we are in range according to the direction we are going
failed = 0;
while ((exp.cont_par.direction >= 0) && (exp.data(end).x_amp < exp.cont_par.x_range(2)) && (exp.data(end).out_amp < exp.cont_par.out_range(2))) || ...
        ((exp.cont_par.direction < 0) && (target > 0) && (exp.data(end).x_amp > exp.cont_par.x_range(1)) && (exp.data(end).out_amp > exp.cont_par.out_range(1)))
    % Set the new target
    target = lasttarget + h;
    rtc.par.x_target_coeffs(idx_sin) = target;
    % Do the Picard iteration
    if ~cbc_picard_iteration_x(exp)
        warning([mfilename ':cbc_picard_iteration_x'], '%s: Picard iteration failed to converge', datestr(now, 13));
    end
    % If adaptive...
    if exp.cont_par.adaptive
        % Check what the step actually was
        x_amp = norm(rtc.par.x_coeffs_ave(exp.fourier.idx_AC));
        out_amp = norm(rtc.par.out_coeffs_ave(exp.fourier.idx_AC));
        step = sqrt((exp.cont_par.step_norm(1)*(exp.data(end).x_amp - x_amp))^2 + (exp.cont_par.step_norm(2)*(exp.data(end).out_amp - out_amp))^2);
        % Apply a linear scaling - probably not reliable for small steps
        h_scale = step_ideal/step;
        if h_scale > exp.cont_par.h_max_scale
            % Arbitrary safety factor
            h_scale = exp.cont_par.h_max_scale;
        end
        h = h*h_scale;
        if (step < exp.cont_par.step_range(1)) || (step > exp.cont_par.step_range(2))
            warning([mfilename ':rejected_step'], '%s: Rejected step - step size is %g', datestr(now, 13), h);
            failed = failed + 1;
            if failed > exp.cont_par.max_fails
                warning([mfilename ':max_rejected_steps'], '%s: Too many rejected steps', datestr(now, 13));
                break
            else
                continue
            end
        end
        failed = 0;
        if exp.cont_par.debug
            fprintf('Step size %g\n', h);
        end
    end        
    lasttarget = target;
    % Record the data
    pt = cbc_get_data_point(exp);
    % Add in the extra data we calculate
    pt.x_amp = norm(pt.x_coeffs_ave(exp.fourier.idx_AC));
    pt.out_amp = norm(pt.out_coeffs_ave(exp.fourier.idx_AC));
    exp.data = [exp.data pt];
    % Save the data
    exp.timestamp_last = clock();
    save(fullfile('data', exp.name), 'exp');
    % Plot the data
    if exp.cont_par.plotting
        plot(ax, exp.data(end).out_amp, exp.data(end).x_amp, [exp.cont_par.colour '.']);
        title(ax, sprintf('forcing = %g, response = %g', exp.data(end).out_amp, exp.data(end).x_amp));
        drawnow;
    end
end

end