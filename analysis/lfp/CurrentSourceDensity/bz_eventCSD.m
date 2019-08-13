
function [ csd, lfpAvg ] = bz_eventCSD (lfp, events, varargin)

% [ CSD ] = bz_eventCSD (lfp, events, varargin)
% Calculates event-triggered (i.e. SWRs) CSD map from a linear array of LFPs

% INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%   events          events timestamps (in sec)

%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       channels    vector with channels to inlcude. If empty take all (default)
%       twin        time window around events to calculate average. Default [0.1 0.1]
%       spat_sm     degree of spatial smoothing. Default = 0.
%       temp_sm     degree of temporal smoothing. Default = 0.
%       plotCSD     true/false. Default true.
%       plotLFP     true/false. Default true.
%    =========================================================================

% OUTPUT:
%    CSD            a buzcode structure with fields csd.data,
%                                                   csd.timestamps
%                                                   csd.samplingRate
%                                                   csd.channels 
%                                                   csd.params
%    lfpAvg         a buzcode structure with fields lpf.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                                                   lfp.channels 
%                                                   lfp.params

% Antonio FR, 7/18

% TODO: make so it can read from the binary lfp file in chunks instead of
% from a lfp.mat
%

%% Parse inputs

p = inputParser;
addParameter(p,'channels',1:size(lfp.data,2),@isvector);
addParameter(p,'samplingRate',1250,@isnumeric);
addParameter(p,'twin',[0.1 0.1],@isnumeric);
addParameter(p,'spat_sm',11,@isnumeric);
addParameter(p,'temp_sm',11,@isnumeric);
addParameter(p,'doDetrend',false,@islogical);
addParameter(p,'plotCSD',true,@islogical);
addParameter(p,'plotLFP',true,@islogical);

parse(p,varargin{:});
channels = p.Results.channels;
samplingRate = p.Results.samplingRate;
spat_sm = p.Results.spat_sm;
temp_sm = p.Results.temp_sm;
doDetrend = p.Results.doDetrend;
plotCSD = p.Results.plotCSD;
plotLFP = p.Results.plotLFP;

%lfp input
if isstruct(lfp)
    data = lfp.data;
    timestamps = lfp.timestamps;
    samplingRate = lfp.samplingRate;
elseif iscell(lfp) %for multiple trials
    celllengths = cellfun(@length,lfp);
    data = vertcat(lfp{:});
elseif isnumeric(lfp)
    data = lfp;
    timestamps = [1:length(lfp)]'./samplingRate;
end

twin = p.Results.twin*samplingRate;
events = round(events*samplingRate);

%% Conpute event-triggered LFP average

events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
lfp_temp = nan(twin(1)+twin(2)+1,length(channels),length(events));

for e = 1:length(events)
    lfp_temp(:,:,e) = data(events(e)-twin(1):events(e)+twin(2),channels);
end

lfp_avg = nanmean(lfp_temp,3)*-1;

%% Conpute CSD

% detrend
if doDetrend
   lfp_avg = detrend(lfp_avg')';
end
    
% temporal smoothing
if temp_sm > 0
   for ch = 1:size(lfp_avg,2) 
       lfp_avg(:,ch) = smooth(lfp_avg(:,ch),temp_sm,'sgolay');
   end
end

% spatial smoothing
if spat_sm > 0
   for t = 1:size(lfp_avg,1) 
       lfp_avg(t,:) = smooth(lfp_avg(t,:),spat_sm,'lowess');
   end
end

% calculate CSD 
CSD = diff(lfp_avg,2,2);

% generate output structure
csd.data = CSD;
csd.timestamps = -twin(1):twin(2);
csd.samplingRate = samplingRate;
csd.channels = channels; 
csd.params.spat_sm = spat_sm;
csd.params.temp_sm = temp_sm;
csd.params.detrend = doDetrend;

lfpAvg.data = lfp_avg;
lfpAvg.timestamps = -twin(1):twin(2);
lfpAvg.samplingRate = samplingRate;
lfpAvg.channels = channels; 
lfpAvg.params.spat_sm = spat_sm;
lfpAvg.params.temp_sm = temp_sm;
lfpAvg.params.detrend = doDetrend;

%% Plot

if plotLFP
    
    taxis = (-(twin(1)/samplingRate):(1/samplingRate):(twin(2)/samplingRate))*1e3;
    cmax = max(max(CSD)); 
    figure;
    subplot(1,2,1);
    contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
    colormap jet; caxis([-cmax cmax]);
    set(gca,'YDir','reverse');xlabel('time (s)');ylabel('channel');title('CSD'); 
    plot([0 0],[1 size(CSD,2)],'--k');hold on;
    
    subplot(1,2,2);
    for ch=1:size(lfp_avg,2)
        offset = 400*(ch-1);
        sh_tmp = 1e0*(lfp_avg(:,ch)) + offset;
        plot(taxis,sh_tmp,'k','LineWidth',1.5); hold on;
        clear sh_tmp
    end
    set(gca,'YDir','reverse','YTickLabel',[]);ylim([-1000 offset+1000]);xlim([taxis(1) taxis(end)]);
    xlabel('time (ms)');ylabel('channel');title('LFP');   
    plot([0 0],ylim,'--r');hold on;

       
elseif plotCSD  
    
     cmax = max(max(CSD)); 
   
     figure;
     contourf(taxis,1:size(CSD,2),CSD',40,'LineColor','none');hold on;
     colormap jet; caxis([-cmax cmax]);
     set(gca,'YDir','reverse','YTickLabel',[]);ylim([-1000 offset+1000]);xlim([taxis(1) taxis(end)]);
     plot([0 0],[1 size(CSD,2)],'--k');hold on;
   
end

end
