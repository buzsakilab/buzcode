function [popBursts] = bz_FindPopBursts(varargin)
% USAGE
%  [popBursts] = bz_FindPopBursts(spikes)
%
% INPUTS
%   
%   spikes
%   threshold       # of stds above mean
%   durations       window for event durations, durations(1)*2 is also ISI minimum
%  
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
%-add detector info and params
%
% this function detects population burst events as described in Skaggs et
% al XXXX
%
% written by david tingley 2018

p = inputParser;
addRequired(p,'spikes',@isstruct)
addParameter(p,'threshold',[3],@isnumeric)
addParameter(p,'speedThresh',[],@isnumeric) % assumes cm/s
addParameter(p,'behavior',[],@isstruct)
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
saveMat = p.Results.saveMat;
speedThresh = p.Results.speedThresh;
behavior = p.Results.behavior;

%%

spkMat = bz_SpktToSpkmat(spikes,'binSize',binSize,'overlap',1);

% save memory and add to single vector...
popRate = zeros(length(spkMat.data(:,1)),1);
for spk = 1:length(spikes.times)
    popRate = popRate + zscore(Smooth(spkMat.data(:,spk),smoothingWindow));
end
popRate = popRate./spk;

[pks locs width] = findpeaks(popRate,'minpeakheight',mean(popRate)+threshold*std(popRate),...
    'WidthReference','halfprom','MinPeakDistance',durations(1)*2);

 
if ~isempty(speedThresh)  % exclude events when speed is above thresh
    exclude = [];
    for pk = 1:length(locs)
       [a b] = min(abs(spkMat.timestamps(locs(pk))-behavior.timestamps));
       if b <=length(behavior.velocity)
       vel = behavior.velocity(b);
       if vel > speedThresh
           exclude = [exclude; pk];
       end
       end
    end
    disp(['excluding ' num2str(length(exclude)) '/' num2str(length(locs)) ' for speed threshold'])
    pks(exclude) = [];
    locs(exclude) = [];
    width(exclude) = [];
end


keep = find(width<durations(2) & width>durations(1));
pks = pks(keep);
locs = locs(keep);
width = width(keep);

starts = nan(length(width),1);
stops = nan(length(width),1);

exclude=[];
for event = 1:length(pks)
    if locs(event)-durations(2) > 0
        back = durations(2) - find(popRate(locs(event)-durations(2):locs(event))<0,1,'last');
    else
        back = durations(2) - find(popRate(1:locs(event))<0,1,'last');
    end
    if locs(event)+durations(2) <= length(popRate)
        forward = find(popRate(locs(event):locs(event)+durations(2))<0,1,'first');
    else
        forward = find(popRate(locs(event):end)<0,1,'first')-1;
    end
   
    if ~isempty(back) & ~isempty(forward) & back < locs(event) 
        starts(event) = locs(event)-back;
        stops(event) = locs(event)+forward;
    else
        exclude = [exclude;event];
        starts(event) = nan;
        stops(event) = nan;
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
popBursts.width = [stops-starts];
popBursts.amplitudes = pks;
popBursts.nSpikes = [];
popBursts.meanSpikes = [];
for event = 1:length(popBursts.amplitudes)
    popBursts.nSpikes(event) = sum(sum(spkMat.data(starts(event):stops(event),:)));
    popBursts.meanSpikes(event) = mean(mean(spkMat.data(starts(event):stops(event),:)));
end
popBursts.nSpikes = popBursts.nSpikes';
popBursts.meanSpikes = popBursts.meanSpikes';



if saveMat
    save([popBursts.sessionName '.popBursts.events.mat'],'popBursts')
end







