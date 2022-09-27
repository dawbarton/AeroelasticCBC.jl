function data = cbc_continex_create_experiment(rtc, varargin)

% Parse the input
p = inputParser();
p.KeepUnmatched = true;
p.addRequired('par', @iscellstr);
p.parse(varargin{:});

% Create a continex specific data structure
continex.par = p.Results.par;

% Go through and extract scalings for different parameters
% - scalings are assumed to be polynomials in a form suitable for polyval
%   or a function handle to do the scaling
unmatched = p.Unmatched;
for par = continex.par
    par_name = par{1};
    field_name = [par_name, '_scale'];
    if isfield(unmatched, field_name)
        % Scaling
        scale = unmatched.(field_name);
        if numel(scale) == 1
            scale = [scale(1), 0];
        elseif numel(scale) > 2
            error('Scaling can only be a linear polynomial');
        end
        inv_scale = [1/scale(1), -scale(2)/scale(1)];
        % Store the scaling as a polyval call
        continex.par_scale.(par_name) = @(x)(polyval(scale, x));
        continex.par_inv_scale.(par_name) = @(x)(polyval(inv_scale, x));
        % Remove the scaling from the unmatched fields
        unmatched = rmfield(unmatched, field_name);
    else
        % Default identity scaling
        continex.par_scale.(par_name) = @(x)x;
        continex.par_inv_scale.(par_name) = @(x)x;
    end
end

% Create the container for the experiment
data = cbc_create_experiment(rtc, unmatched);
data.continex = continex;

end
