function [M, x_0, x_1] = poincare_map(rtc, perturb, n_samples, n_cycles_per_sample)

% Size and duration of perturbations
perturb_default = [1, 0.25]; 
if ~exist('perturb', 'var')
    perturb_duration = perturb_default(1);
    perturb_size = perturb_default(2);
else
    if numel(perturb) == 1
        perturb_duration = perturb;
        perturb_size = perturb_default(2);
    else
        perturb_duration = perturb(1);
        perturb_size = perturb(2);
    end
end

% Number of samples (perturbations) to make
if ~exist('n_samples', 'var')
    n_samples = 40;
end

% Number of cycles of decay to consider; for high damping this should be
% low as the later cycles don't contain any information but for low damping
% this can be high(er) since the later cycles are still decaying
if ~exist('n_cycles', 'var')
    n_cycles_per_sample = 2;
end

% Sample frequency for local reference
sample_freq = double(rtc.par.sample_freq);

% Duration for stream recording - get at least n_cycles_per_sample+1 cycles
duration = ceil((n_cycles_per_sample + 1)*double(sample_freq)/rtc.par.forcing_freq) + 1;

% Set the RTC to record what we are interested in
rtc.set_stream(1, {'time_mod_2pi', 'x'}, duration, 0);

% Get the basic motion
data = rtc.run_stream(1, 'struct', true);

% Smooth the data
[B, A] = butter(8, 0.05, 'low');
data.x = filtfilt(B, A, double(data.x));

% Calculate the derivative
data.dx = conv(data.x, [1, 0, -1], 'same')*0.5*sample_freq;
data.dx(1) = (data.x(2) - data.x(1))*sample_freq;
data.dx(end) = (data.x(end) - data.x(end - 1))*sample_freq;

% Fix the time index
data.time_mod_2pi = data.time_mod_2pi/(2*pi) + cumsum([0, double(diff(data.time_mod_2pi) < 0)]);

% Work out where the fixed point is
x_vals = interp1(data.time_mod_2pi', [data.x; data.dx]', 1:n_cycles_per_sample+1);
x_fixed = mean(x_vals)';

% Criteria for being close to the orbit
eps = 3*norm(sqrt(rtc.par.x_coeffs_var));

% Perturb the orbit and measure the response
x_0 = zeros(2, n_samples*n_cycles_per_sample);
x_1 = zeros(2, n_samples*n_cycles_per_sample);
idx = 1;
forcing_amp = rtc.par.forcing_amp;
if forcing_amp == 0
    forcing_amp_perturb = perturb_size;
else
    forcing_amp_perturb = forcing_amp*perturb_size;
end
% forcing_amp = rtc.par.forcing_amp;
for i = 1:n_samples
    % Perturb
    if i <= n_samples/2
        rtc.par.forcing_amp = forcing_amp - forcing_amp_perturb;
        pause(perturb_duration);
        rtc.par.forcing_amp = forcing_amp;
    else
        rtc.par.forcing_coeffs(rtc.fourier.idx_cos(1)) = perturb_size;
        pause(perturb_duration);
        rtc.par.forcing_coeffs(rtc.fourier.idx_cos(1)) = 0;
    end
    % Get the data
    data = rtc.run_stream(1, 'struct', true);
    % Smooth the data
    [B, A] = butter(8, 0.05, 'low');
    data.x = filtfilt(B, A, double(data.x));
    % Calculate the derivative
    data.dx = conv(data.x, [1, 0, -1], 'same')*0.5*sample_freq;
    data.dx(1) = (data.x(2) - data.x(1))*sample_freq;
    data.dx(end) = (data.x(end) - data.x(end - 1))*sample_freq;
    % Fix the time index
    data.time_mod_2pi = data.time_mod_2pi/(2*pi) + cumsum([0, double(diff(data.time_mod_2pi) < 0)]);
    % Find the x values at t=0 (mod 2 pi)
    x_vals = interp1(data.time_mod_2pi', [data.x; data.dx]', 1:n_cycles_per_sample+1)';
    j_max = min([length(x_vals), find(abs(x_vals) > eps, 1, 'last')]);
    % Iterate over the found (and acceptably far away from steady-state) values
    for j = 2:j_max
        x_0(:, idx) = x_vals(:, j - 1) - x_fixed;
        x_1(:, idx) = x_vals(:, j) - x_fixed;
        idx = idx + 1;
    end
end
x_0 = x_0(:, 1:idx-1);
x_1 = x_1(:, 1:idx-1);

figure; plot3(x_0(1, :), x_0(2, :), x_1(1,:), '.');

% Use least-squares to calculate the monodromy matrix
M = x_1 / x_0;

end