function [popBursts] = bz_FindPopBursts(varargin)
% USAGE
%  [popBursts] = bz_FindPopBursts(spikes)
%
% INPUTS
%   
%   spikes
%   threshold       # of stds above mean
%   durations       window for event durations
%
%
% OUTPUTS
%
%
%
%
%
% TODO
% -zscore option
% -param input parser
% - smoothing binSize
% -% get start/stop times 
%
%
% this function detects population burst events as described in Skaggs et
% al XXXX
%
% written by david tingley 2018

p = inputParser;
addRequired(p,'spikes',@isstruct)
addParameter(p,'threshold',[3],@isnumeric)
addParameter(p,'durations',[10 500],@isnumeric)
addParameter(p,'binSize',.001,@isnumeric)
addParameter(p,'smoothingWindow',10,@isnumeric)
addParameter(p,'saveMat',false,@islogical)


parse(p,varargin{:})

spikes = p.Results.spikes;
threshold = p.Results.threshold;
durations = p.Results.durations;
binSize = p.Results.binSize;
smoothingWindow = p.Results.smoothingWindow;

%%

spkMat = bz_SpktToSpkmat(spikes,'binSize',binSize,'overlap',1);

for spk = 1:length(spikes.times)
    spkMat.zscoredData(:,spk) = zscore(Smooth(spkMat.data(:,spk),smoothingWindow));
end

popRate = nanmean(spkMat.zscoredData');

[pks locs width] = findpeaks(popRate,'minpeakheight',mean(popRate)+threshold*std(popRate),...
    'WidthReference','halfprom');

keep = find(width<durations(2) & width>durations(1));
pks = pks(keep);
locs = locs(keep);
width = width(keep);

exclude=[];
for event = 1:length(pks)
    forward = find(popRate(locs(event):locs(event)+durations(2))<0,1,'first');
    back = durations(2) - find(popRate(locs(event)-durations(2):locs(event))<0,1,'last');
    if ~isempty(back) & ~isempty(forward)
        starts(event) = locs(event)-back;
        stops(event) = locs(event)+forward;
    else
        exclude = [exclude;event];
    end
end

pks(exclude) = [];
locs(exclude) = [];
width(exclude) = [];
starts(exclude) = [];
stops(exclude) = [];

popBursts.sessionName = spikes.sessionName;
popBursts.timestamps = [spkMat.timestamps(starts) spkMat.timestamps(stops)];
popBursts.bursts = spkMat.timestamps(locs);
popBursts.width = [stops-starts]';




if saveMat
    save([popBursts.sessionName '.popBursts.events.mat'],'popBursts')
end
c













