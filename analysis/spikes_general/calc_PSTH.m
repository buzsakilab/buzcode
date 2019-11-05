function PSTH = calc_PSTH(event,spikes,varargin)
% This is a generalized way for creating a PSTH for units for various events
% 
% INPUTS
% event  : event times formatted according to the Cell Explorer's convention
% spikes : spikes formatted according to the Cell Explorer's convention
% 
% OUTPUT
% psth
% 
% Dependencies: CCG

% By Peter Petersen
% petersen.peter@gmail.com
% Last edited 11-08-2019

p = inputParser;

addParameter(p,'binCount',200,@isnumeric);        % how many bins (for half the window)
addParameter(p,'alignment','onset',@ischar);    % alignment of time ['onset','center','peaks','offset']
addParameter(p,'binDistribution',[0.25,0.5,0.25],@isnumeric);  % How the bins should be distributed around the events, pre, during, post. Must sum to 1
addParameter(p,'duration',0,@isnumeric);        % duration of PSTH (for half the window - used in CCG) [in seconds]
addParameter(p,'smoothing',0,@isnumeric);       % any gaussian smoothing to apply? units of bins.
addParameter(p,'percentile',99,@isnumeric);     % if events does not have the same length, the event duration can be determined from percentile of the distribution of events
addParameter(p,'eventName','',@ischar); 
addParameter(p,'plots',true,@islogical);

parse(p,varargin{:})

binCount = p.Results.binCount;
alignment = p.Results.alignment;
binDistribution = p.Results.binDistribution;
duration = p.Results.duration;
smoothing = p.Results.smoothing;
percentile = p.Results.percentile;
eventName = p.Results.eventName;
plots = p.Results.plots;

% If no duration is given, an optimal duration is determined
if duration == 0
    durations = diff(event.timestamps');
    stim_duration = prctile(sort(durations),percentile);
    duration = min(max(round(stim_duration*1000),50)/1000,30);
end

binSize = max(round(duration/binCount*1000),1)/1000; % minimum binsize is 0.5ms.

% Determine event alignment
switch alignment
    case 'onset'
        event_times = event.timestamps(:,1);
        padding = binDistribution(1)/binDistribution(2)*stim_duration;
        binsToKeep = ceil((padding+duration/2)/binSize):(duration+padding)*2/binSize;
    case 'center'
        event_times = mean(event.timestamps);
        padding = 0;
        binsToKeep = 1:duration*2/binSize;
    case 'offset'
        event_times = event.timestamps(:,2);
        padding = binDistribution(3)/binDistribution(2)*stim_duration;
        binsToKeep = 1:(duration+padding)*2/binSize-ceil((padding+duration/2)/binSize);
    case 'peaks'
        event_times = event.peaks;
        padding = 0;
        binsToKeep = 1:duration*2/binSize;
end

disp(['  ', num2str(length(event_times)), '  events, duration set to: ', num2str(duration), ' sec, aligned to ', alignment])

% Determining the bins interval for metrics
binsPre = 1:floor(binDistribution(1)*length(binsToKeep));
binsEvents = floor(binDistribution(1)*length(binsToKeep))+1:floor((binDistribution(1)+binDistribution(2))*length(binsToKeep));
binsPost = floor((binDistribution(1)+binDistribution(2))*length(binsToKeep))+1:length(binsToKeep);

% Calculating PSTH
spike_times = spikes.spindices(:,1);
spike_cluster_index = spikes.spindices(:,2);
[spike_times,index] = sort([spike_times;event_times(:)]);
spike_cluster_index = [spike_cluster_index;zeros(length(event_times),1)];
spike_cluster_index = spike_cluster_index(index);
[~, ~, spike_cluster_index] = unique(spike_cluster_index);
[ccg,time] = CCG(spike_times,spike_cluster_index,'binSize',binSize,'duration',(duration+padding)*2);

time = time(binsToKeep+1);
PSTH_out = flip(ccg(:,2:end,1),1);
% PSTH_out = ccg(:,2:end,1);
PSTH_out = PSTH_out(binsToKeep+1,:)./length(event_times)/binSize;

modulationIndex = mean(PSTH_out(binsEvents,:))./mean(PSTH_out(binsPre,:));
modulationSignificanceLevel = [];
for i = 1:size(PSTH_out,2)
    [~,p_kstest2] = kstest2(PSTH_out(binsEvents,i),PSTH_out(binsPre,i));
    modulationSignificanceLevel(i) = p_kstest2;
end

if smoothing>0
    PSTH_out = nanconv(PSTH_out,gausswin(smoothing)/sum(gausswin(smoothing)),'edge');
end

[~,modulationPeakResponseTime] = max(PSTH_out);
modulationPeakResponseTime = time(modulationPeakResponseTime);

PSTH.responsecurve = PSTH_out;
PSTH.time = time;
PSTH.alignment = alignment;

PSTH.modulationIndex = modulationIndex;
PSTH.modulationPeakResponseTime = modulationPeakResponseTime';
PSTH.modulationSignificanceLevel = modulationSignificanceLevel;

if plots
    figure, plot(time,PSTH_out), title(eventName), xlabel('Time')
    [~,index2] = sort(modulationIndex,'descend');
    [~,index3] = sort(modulationPeakResponseTime);
    
    figure,
    subplot(2,2,1), histogram(modulationIndex,40), title('modulationIndex'), xlabel('Ratio'), ylabel(eventName)
    subplot(2,2,2), histogram(modulationPeakResponseTime,40), title('modulationPeakResponseTime'), xlabel('Time')
    subplot(2,2,3), imagesc(time,[1:size(PSTH_out,2)],zscore(PSTH_out(:,index2))'), title(['Sorting: modulationIndex']), xlabel('Time'), ylabel('Units')
    subplot(2,2,4), imagesc(time,[1:size(PSTH_out,2)],zscore(PSTH_out(:,index3))'),  title(['Sorting: modulationPeakResponseTime']), xlabel('Time')
end
