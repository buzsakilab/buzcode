function [spectrum, f] = long_wavespec(signal, t_dt_sec, window_sec, frange, df)

if (length(t_dt_sec) > 1)
    dt = diff(t_dt_sec(1:2));
    t_sec = t_dt_sec;
else
    dt = t_dt_sec;
    t_sec = (0:dt:((length(signal)-1)*dt))';
end
Fs = 1/dt;

nWind = window_sec/dt + 1;
wind_overlap = 2*round(0.05*nWind/2); % force even

% Set window intervals
% First window only overlaps on one side; accounting for this, find the
% number of windows based on the length of the signal
whole_winds = 1 + floor((length(t_sec) - nWind + wind_overlap/2)/(nWind - wind_overlap));
winds = cell(whole_winds + 1, 1);
winds{1} = 1:nWind;
for i=2:length(winds)
    winds{i} = winds{i-1}(end) - wind_overlap + (1:nWind);
end
winds{end} = winds{end}(winds{end} <= length(t_sec));
% If final window shorter than three periods of lowest freqency of
% interest, do not include it in the analysis
if ((length(winds{end}) - 1)*dt < 3/frange(1))
    winds = winds(1:(end-1));
end

% frequencies
f = frange(1):df:frange(2);
cwtScales = Fs./f;

spectrum = zeros(length(f), length(signal));
for i=1:length(winds)
    complexCWT_wind = cwt(signal(winds{i}), cwtScales, 'cmor1.5-1');
    spectrum(:,winds{i}((1 + wind_overlap/2):(end - wind_overlap/2))) = abs(complexCWT_wind(:,(1 + wind_overlap/2):(end - wind_overlap/2)));
end
