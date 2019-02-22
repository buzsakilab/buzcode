function [cell_metrics] = bz_CellMetricsSimple(spikes,sr)
% [cell_metrics] = bz_CellMetricsSimple(spikes,sr)
%   Calculate waveform and ACG metrics to use for unit classification
%
% USAGE
%
%   [cell_metrics] = bz_CellMetricsSimple(spikes)
%
% INPUTS
%   
%   spikes = buzcode spikes structure (i.e. from bz_LoadPhy)
%   sr = sampling rate
%
% OUTPUT
%
%   cell_metrics
%       .filtWaveform
%       .rawWaveform
%           .PeaktoTrough 
%           .TroughtoPeak
%           .AB_ratio
%           .trough
%       .refractory
%       .burstRoyer
%       .burstiness
%       .FR

% Antonio FR, 11/2018

if ~isnumeric(sr)
    sr = 30000;
end

%% waveform metrics
   % reshape data 
   for i = 1:length(spikes.rawWaveform)
       rawWaveforms(i,:) = spikes.rawWaveform{i};
       filtWaveforms(i,:) = spikes.filtWaveform{i};
   end
   
   wave = [];t_before = [];t_after = [];peakA = [];peakB = [];trough = [];
   for m = 1:size(rawWaveforms,1)
        wave = interp1([1:size(filtWaveforms,2)],zscore(filtWaveforms(m,:)),[1:0.5:size(filtWaveforms,2),size(filtWaveforms,2)],'spline');
        midPoint = round((length(wave)/2));
        [MIN2,I2] = min(wave(midPoint-10:midPoint+10));
        [MAX3,I3] = max(wave(1:midPoint));
        [MAX4,I4] = max(wave((I2+midPoint-11):end));
        t_before(m) = (I2+midPoint-11)-I3;
        t_after(m) = I4;
        peakA(m) = MAX3;
        peakB(m) = MAX4;
        trough(m) = MIN2;
        clear MIN2 MAX3 MAX4 I2 I3 I4
   end
    cell_metrics.filtWaveform.PeaktoTrough = (t_before/(sr*2))';
    cell_metrics.filtWaveform.TroughtoPeak = (t_after/(sr*2))';
    cell_metrics.filtWaveform.AB_ratio = ((peakB-peakA)./(peakA+peakB))';
    cell_metrics.filtWaveform.trough = (trough)';      
   
   wave = [];t_before = [];t_after = [];peakA = [];peakB = [];trough = [];
   for m = 1:size(rawWaveforms,1)
        wave = interp1([1:size(rawWaveforms,2)],zscore(rawWaveforms(m,:)),[1:0.5:size(rawWaveforms,2),size(rawWaveforms,2)],'spline');
        midPoint = round((length(wave)/2));
        [MIN2,I2] = min(wave(midPoint-10:midPoint+10));
        [MAX3,I3] = max(wave(1:midPoint));
        [MAX4,I4] = max(wave((I2+midPoint-11):end));
        t_before(m) = (I2+midPoint-11)-I3;
        t_after(m) = I4;
        peakA(m) = MAX3;
        peakB(m) = MAX4;
        trough(m) = MIN2;
        clear MIN2 MAX3 MAX4 I2 I3 I4
   end
    cell_metrics.rawWaveform.PeaktoTrough = (t_before/(sr*2))';
    cell_metrics.rawWaveform.TroughtoPeak = (t_after/(sr*2))';
    cell_metrics.rawWaveform.AB_ratio = ((peakB-peakA)./(peakA+peakB))';
    cell_metrics.rawWaveform.trough = (trough)';   
     
%%  ACG metrics     
for i = 1:length(spikes.times)
    [ccg,time] = CCG(spikes.times{i},ones(length(spikes.times{i}),1),'binSize',0.0005,'duration',0.100); %100ms wide CCG with 0.5ms bins
    ccg0 = find(time == 0);ccg10 = find(time == 0.01);ccg20 = find(time == 0.02);
    ccg40 = find(time == 0.04);ccg50 = find(time == 0.05);

    % Refractory period and burst index Royer 2012
    CCGmax = find(ccg(ccg0:ccg20)== max(ccg(ccg0:ccg20)));
    CCGmax = CCGmax+ccg0-1;  
    if ~isempty(find(diff(ccg(ccg0:CCGmax(1)))>std(diff(ccg(ccg0:CCGmax(1)))),1))% refractory period according to Royer 2012
        cell_metrics.refractory(i,1) = find(diff(ccg(ccg0:CCGmax(1)))>std(diff(ccg(ccg0:CCGmax(1)))),1);
    else
        cell_metrics.refractory(i,1) = NaN;
    end

    CCGbase=mean(ccg(ccg40:ccg50));clear CCGmax; % mean CCG 40-50 ms
    CCGmax=find(ccg(ccg0:ccg10)==max(ccg(ccg0:ccg10)));
    CCGmax = CCGmax+ccg0-1; % only from 0 to 10ms

    if ccg(CCGmax)>CCGbase
       cell_metrics.burstRoyer(i,1)=(ccg(CCGmax(1))-CCGbase)/ccg(CCGmax(1));
    else
       cell_metrics.burstRoyer(i,1)=(ccg(CCGmax(1))-CCGbase)/(CCGbase); 
    end
    
    % Burstiness 
     binssum = 6; %for isis below 6ms 
     binSize = 1e-3; nBins = 100; binEdges = (0:nBins) * binSize;
     t = ((0:nBins-1) + 0.5)' * binSize;
     ISIs = diff(spikes.times{i});
     hISI = histc(ISIs, binEdges); 
     cell_metrics.burstiness(i,1) = sum(hISI(1:binssum))/sum(hISI);   

    % Sam's ACG fit 
    rsquare = [];   offset = 106;
    g = fittype('c*exp(-x/a)-d*exp(-x/b)');
    for j = 1:size(ccg,2)
        x = [1:95]';
        y = ccg(x+offset,j)/max(ccg(x+offset,j));

        [f0,gof] = fit(x/2,y,g,'StartPoint',[25, 1, 1, 1.5],'Lower',[1,0.1,-Inf, -Inf],'Upper',[100, 10,2, 10]);
    %     f0 = fit(x/2,y,'exp2','StartPoint',[1, -0.015, -1, -1]);
        fit_params(:,j) = coeffvalues(f0);
        xx = linspace(0,48,100);
        rsquare(j) = gof.rsquare;
        %figure, plot(x/2,y,'-o',x/2,f0(x/2),'r-'); title(['Unit ', num2str(j), ' r^2: ' num2str(rsquare(j))])
    end

    ACG.tau_decay = fit_params(1,:);
    ACG.tau_rise = fit_params(2,:);
    ACG.c = fit_params(3,:);
    ACG.d = fit_params(4,:);
    ACG.fit_rsquare = rsquare;
    
    clear s ccg CCGmax refT CCGbase ISIs hISI ; 
end

%% FR
for i = 1:length(spikes.times)
    cell_metrics.FR(i,1) = length(spikes.times{i})/(spikes.times{i}(end)-spikes.times{i}(1));
end

%% Preliminary pyr/int classification
    % avg FR > 6Hz or narrow waveform (trought to peak < 5 ms) = INT
for i = 1:length(spikes.times)
    if cell_metrics.FR(i,1) > 6 || cell_metrics.filtWaveform.TroughtoPeak(i) < 5e-4
        cell_metrics.putativeClass(i,1) = 2; % int
    else 
        cell_metrics.putativeClass(i,1) = 1; % pyr
    end
end

end

%% 
