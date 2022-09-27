function exp = cbc_create_experiment(rtc, varargin)
% CBC_CREATE_EXPERIMENT  Create a data structure intended to hold the
%   output of a particular experiment run.
%
%   EXP = CBC_CREATE_EXPERIMENT(RTC) returns an experiment data structure,
%   based on the interface to the real-time controller RTC, that can be
%   used in other CBC functions, such as CBC_CONTINUATION_AMPLITUDE. 
%
%   EXP = CBC_CREATE_EXPERIMENT(RTC, Name, Value) overrides any of the
%   default options.
%
%   Static data parameters (indicated by the static_fields option) are
%   captured when the experiment is created. Dynamic and stream data
%   parameters (indicated by the dynamic_fields and stream_fields options
%   respectively) are captured on subsequent calls to CBC_GET_DATA_POINT.
%
%   By default the static_fields, dynamic_fields and stream_fields options
%   add fields to be captured in addition to the defaults. To remove the
%   defaults set the option default_fields to false.
%
%   Any options in the RTC.opts field can also be overridden.
%
%   See also CBC_CONTINUATION_AMPLITUDE and CBC_GET_DATA_POINT.

% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015
 
% Parse the input
p = inputParser();
if ismethod(p, 'addParameter')
    % New versions of Matlab
    add_par = @p.addParameter;
else
    % Old versions of Matlab
    add_par = @p.addParamValue; %#ok<NVREPL>
end
add_par('default_fields', true, @islogical);
add_par('static_fields', {}, @iscellstr);
add_par('dynamic_fields', {}, @iscellstr);
add_par('stream_fields', {}, @iscellstr);
cellfun(add_par, fieldnames(rtc.opt), struct2cell(rtc.opt)); % Add all fields in opt as potential options
p.parse(varargin{:});

% Store the interface to the RTC
exp.rtc = rtc;

% Name of the experiment (to be set by the user); typically a file name
exp.name = '';

% Time stamp when the experiment was created/last modified
exp.timestamp_created = clock();
exp.timestamp_last = exp.timestamp_created;

% Record the settings of the experiment
exp.fourier = rtc.fourier;
exp.opt = rtc.opt;
exp.datafields = rtc.datafields;

% Options for data recording
exp.opt.samples = p.Results.samples; % Number of samples to take
exp.opt.downsample = p.Results.downsample; % Number of samples to ignore for each sample taken

% Default fields to record
if p.Results.default_fields
    exp.datafields.static_fields = [exp.datafields.static_fields, p.Results.static_fields];
    exp.datafields.dynamic_fields = [exp.datafields.dynamic_fields, p.Results.dynamic_fields];
    exp.datafields.stream_fields = [exp.datafields.stream_fields, p.Results.stream_fields];
else
    exp.datafields.static_fields = p.Results.static_fields;
    exp.datafields.dynamic_fields = p.Results.dynamic_fields;
    exp.datafields.stream_fields = p.Results.stream_fields;
end

% Record (hopefully!) static parameters
for fieldcell = exp.datafields.static_fields
    field = fieldcell{1};
    exp.par.(field) = rtc.par.(field);
end

% Create the basic data structure for recorded data points
exp.data = [];

% Set RTC stream 1 to record the desired variables
exp.rtc.set_stream(exp.datafields.stream_id, exp.datafields.stream_fields, exp.opt.samples, exp.opt.downsample);

end
