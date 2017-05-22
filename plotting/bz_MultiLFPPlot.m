function [  ] = bz_MultiLFPPlot( lfp,varargin )
%bz_MultiLFPPlot(lfp). This function plots multiple lfp channels from a 
%buzcode lfp structure with the channels appropriately spaced. 
%
%USAGE
%   figure
%       bz_MultiLFPPlot(lfp)
%
%INPUT
%   lfp         a buzcode lfp structure. channels are assumed to be ordered
%               from top to bottom
%   (optional)
%   'channels'  subset/ordering of channels to plot (0-index a la neuroscope)
%               default is to take the channels as ordered in lfp.channels
%   'timewin'   only plot a subwindow of time
%
%
%DLevenstein 2017
%% parse the inputs!
channelsValidation = @(x) assert(isnumeric(x) || strcmp(x,'all'),...
    'channels must be numeric or "all"');

% parse args
p = inputParser;
addParameter(p,'channels','all',@isnumeric)
addParameter(p,'timewin',[0 Inf],@isnumeric)
parse(p,varargin{:})
timewin = p.Results.timewin;
channels = p.Results.channels;

%% Channel and time stuff
%Time Window
windex = lfp.timestamps>=timewin(1) & lfp.timestamps<=timewin(2);
%Channel to data array index mapping
if strcmp(channels,'all')
    chindex = 1:length(lfp.channels);
    channels = lfp.channels;
else
    [~,~,chindex] = intersect(lfp.channels,channels,'stable');
end

%% Calculate and implement spacing between channels

%Space based on median absolute deviation - robust to outliers.
channelrange = 6.*mad(single(lfp.data(windex,chindex)),1);
lfpmidpoints = -cumsum(channelrange);
lfpplotdata = int16(bsxfun(@(X,Y) X+Y,single(lfp.data(windex,chindex)),lfpmidpoints));

%% Do the plot
plot(lfp.timestamps(windex),lfpplotdata,'k')
xlabel('t (s)')
ylabel('LFP by Depth')
set(gca,'Ytick',fliplr(lfpmidpoints))
set(gca,'yticklabels',fliplr(channels))
ylim(fliplr(lfpmidpoints([1 end])+1.*[1 -1].*max(channelrange)))


end

