function pt = cbc_get_data_point(exp)
% CBC_GET_DATA_POINT  Record a data point.
%
%   PT = CBC_GET_DATA_POINT(EXP) measures and returns a data point using
%   the experiment options specified in EXP.
%
%   See also CBC_CREATE_EXPERIMENT.

% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015

% Get the handle to the RTC
rtc = exp.rtc;

% Start the stream running while we are getting the dynamic data
rtc.start_stream(exp.datafields.stream_id);

% Get dynamic data
pt = [];
for fieldcell = exp.datafields.dynamic_fields
    field = fieldcell{1};
    pt.(field) = rtc.par.(field);
end

% Get stream data 
data = rtc.run_stream(exp.datafields.stream_id, 'start', false, 'struct', true);
for fieldcell = exp.datafields.stream_fields
    field = fieldcell{1};
    pt.(field) = data.(field);
end

% Time stamp the data
pt.timestamp = clock();

end
