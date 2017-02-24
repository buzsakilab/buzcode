function [n, npw, ts, wave] = plx_waves(filename, channel, unit)
% plx_waves(filename, channel, unit): read waveform data from a .plx or .pl2 file
%
% [n, npw, ts, wave] = plx_waves(filename, channel, unit)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based channel number or channel name
%   unit  - unit number (0- unsorted, 1-4 units a-d)
%
% OUTPUT:
%   n - number of waveforms
%   npw - number of points in each waveform
%   ts - array of timestamps (in seconds) 
%   wave - array of waveforms [npw, n], raw a/d values

n = 0;
npw = 32;
ts = -1;
wave = -1;

if nargin ~= 3
    error 'expected 3 input arguments';
end

[ filename, isPl2 ] = internalPL2ResolveFilenamePlx( filename );
if isPl2 == 1
    w = PL2Waves( filename, channel, unit );
    if numel(w.Ts) > 0
        n = numel(w.Ts);
        ts = w.Ts;
        npw = w.NumPointsWave;
        wave = w.Waves;
        % convert millivolts to raw a/d values
        pl2 = PL2GetFileIndex(filename);
        channelNumber = internalPL2ResolveChannel(pl2.SpikeChannels, channel);
        wave = round( wave / (pl2.SpikeChannels{channelNumber}.CoeffToConvertToUnits * 1000) );  
    else
		% we need to return number of points per wave the way it is done by mexPlex
		if channel > -1
			npw = 0;
			pl2 = PL2GetFileIndex(filename);
			% find the first channel with data and get its number of points per wave
			for i=1:numel(pl2.SpikeChannels)
				if pl2.SpikeChannels{i}.SamplesPerSpike > 0
					npw = pl2.SpikeChannels{i}.SamplesPerSpike;
					break;
				end
			end
		end
    end
    return;
end

channelNumber = plx_resolve_channel(filename, channel);
if channelNumber == -1
    %fprintf('\n plx_waves: no header for the specified spike channel.\n');
    return
end

[n, npw, ts, wave] = mexPlex(6, filename, channelNumber, unit);

end
