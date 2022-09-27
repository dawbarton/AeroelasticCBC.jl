function [data, y] = cbc_continex_wrapper(data, x, p)

% Set the control target and the desired parameter value
data.rtc.par.x_target_coeffs(data.fourier.idx_fund) = x;
if iscellstr(data.continex.par)
    p = num2cell(p);
    for i = 1:length(data.continex.par)
        p{i} = data.continex.par_inv_scale.(data.continex.par{i})(p{i});
    end
    data.rtc.set_par(data.continex.par, p);
end

% Do a fixed-point (picard) iteration to eliminate the higher modes
cbc_picard_iteration_x(data);

% Get the steady-state value
y = x(:) - data.rtc.par.x_coeffs_ave(data.fourier.idx_fund)';

end
