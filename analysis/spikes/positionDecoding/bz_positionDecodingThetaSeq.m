function [] = thetaSeqDecoder(spikes,lfp,behavior)
% USAGE
%
%
% INPUTS
% 
%   spikes          cell array of spike times where each element is a neuron and
%                   each timestamp is represented in seconds
%   positions       Nx2 matrix of X and Y positions
%   lfp             
%   behavior.events.trialIntervals
%
%
% OUTPUTS
%
%
%
%
%
% HELP
% This function tries to reproduce the theta sequence detection and
% decoding done here: Wikenheiser and Redish. Hippocampal
% theta sequences reflect current goals 2015
% 
% written by david tingley, 2017

numSpikesThresh = 2;

%% first let's get the ratemaps for these cells, and find the place field 'centers'

disp('getting ratemaps and place field statistics..')
if strcmp(behavior.trackingType,'optitrack')
    tau = 5;
    maxFieldWidth = 100;
    minFieldWidth = 8;
else
   tau = 2; 
   maxFieldWidth = 40;
   minFieldWidth = 3;
end
[rateMap countMap opk_rate_vel_corruMap phaseMap] = bz_firingMap1D(spikes.times,behavior,lfp,tau);
    
% remove high firing rate neurons
for i=1:length(spikes.times)
    FR = length(spikes.times{i}) ./ lfp.duration;
    if FR > 5
        spikes.times{i} = [];
    end
end

% let's calculate the population rate
pop_rate = zeros(length(lfp.data),1);
for i=1:length(spikes.times)
    pop_rate(round(1250*spikes.times{i})) = pop_rate(round(1250*spikes.times{i})) + 1;
end

% find place fields and their COM's
for i=1:length(unique(behavior.events.trialConditions))
    fields{i} = bz_getPlaceFields1D(rateMap{i},'minPeakRate',3,'maxFieldWidth',maxFieldWidth,'minFieldWidth',minFieldWidth);
end

%% now let's find all theta thetaCycleQuartilesTS and split them into quartiles
% [b a] = butter(4,[1/(lfp.samplingRate/2) 4/(lfp.samplingRate/2)],'bandpass');
% delta= FiltFiltM(b,a,double(lfp.data(:,1)));
% delta_pow = fastrms(delta,1000);

[b a] = butter(3,[6/(lfp.samplingRate/2) 10/(lfp.samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
pow = fastrms(filt,250);

% phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data(:,1)))));
[blah peaks] = findpeaks(filt);
[blah troughs] = findpeaks(-filt);
% [blah peaks] = findpeaks(cos(phases));
% [blah troughs] = findpeaks(-cos(phases));

disp('segmenting theta cycles into quartiles by peaks/troughs..')
if peaks(1) < troughs(1) % make sure we start on a full cycle (peak to peak)
        peaks(1) = [];
end
    
for i=1:length(peaks) - 1  % parse into quartiles
    thetaCycleQuartilesBins(i,:) = [peaks(i) (peaks(i)+troughs(i+1))/2 troughs(i+1)  (peaks(i+1)+troughs(i+1))/2 peaks(i+1)];
end

% check that there are enough spikes for last theta cycle quartile
keep=[];exclude=[];

for i = 1:size(thetaCycleQuartilesBins,1)
   if sum(pop_rate(floor(thetaCycleQuartilesBins(i,4)):ceil(thetaCycleQuartilesBins(i,5)))) > numSpikesThresh
       keep = [keep; i];
   else
       exclude = [exclude;i];
   end
end

thetaCycleQuartilesBins = thetaCycleQuartilesBins(keep,:);
disp(['searching through ' num2str(size(thetaCycleQuartilesBins,1)) ' candidate theta cycles'])

%% now we calculate the distance between the animals position and the COM 
%% for cells spiking in the last quartile of each theta sequence
thetaCycleQuartilesTS = thetaCycleQuartilesBins ./  lfp.samplingRate;

for t = 1:size(behavior.events.trialIntervals,1)
   trial_thetaCycleQuartilesTS = find(InIntervals(thetaCycleQuartilesTS(:,3),behavior.events.trialIntervals(t,1:2))); 
   if ~isempty(trial_thetaCycleQuartilesTS)
   for tt = 1:length(trial_thetaCycleQuartilesTS)
    	cellIDs = (spikes.spindices((InIntervals(spikes.spindices(:,1),thetaCycleQuartilesTS(trial_thetaCycleQuartilesTS(tt),[4 5]))),2));
        if length(cellIDs) > numSpikesThresh
        for cell = 1:length(cellIDs)
            if ~isempty(fields{behavior.events.trialConditions(t)}{cellIDs(cell)}) & length(fields{behavior.events.trialConditions(t)}{cellIDs(cell)}) == 1
                com(cell) = fields{behavior.events.trialConditions(t)}{cellIDs(cell)}{1}.COM;
                center(cell) = (fields{behavior.events.trialConditions(t)}{cellIDs(cell)}{1}.start+...
                    fields{behavior.events.trialConditions(t)}{cellIDs(cell)}{1}.stop)./2;
            else
                com(cell) = nan;
                
                center(cell) =nan;
            end
        end
%         com(abs(com-1)<10)=nan;
%         com(abs(center-length(behavior.events.map{behavior.events.trialConditions(t)}.x))<10)=nan;
%         center(abs(center-1)<10)=nan;
%         center(abs(center-length(behavior.events.map{behavior.events.trialConditions(t)}.x))<10)=nan;
        predicted_lastquart(tt) = nanmean(center);
        
        [a b] = min(abs(behavior.events.trials{t}.timestamps - thetaCycleQuartilesTS(trial_thetaCycleQuartilesTS(tt),3)));
        actual_pos(tt) = behavior.events.trials{t}.mapping(b);
        vel = smooth(abs(diff(behavior.events.trials{t}.x))+abs(diff(behavior.events.trials{t}.y)),...
            length(behavior.events.map{behavior.events.trialConditions(t)}.x)/10);
        vel = [vel; vel(end)];
        actual_vel(tt) = vel(b);
        clear com center;
        else
            actual_pos(tt)=nan;
            predicted_lastquart(tt) = nan;
            actual_vel(tt) = nan;
        end
        
        cellIDs = (spikes.spindices((InIntervals(spikes.spindices(:,1),thetaCycleQuartilesTS(trial_thetaCycleQuartilesTS(tt),[1 3]))),2));
        if length(cellIDs) > numSpikesThresh
        for cell = 1:length(cellIDs)
            if ~isempty(fields{behavior.events.trialConditions(t)}{cellIDs(cell)}) & length(fields{behavior.events.trialConditions(t)}{cellIDs(cell)}) == 1
                com(cell) = fields{behavior.events.trialConditions(t)}{cellIDs(cell)}{1}.COM;
                center(cell) = (fields{behavior.events.trialConditions(t)}{cellIDs(cell)}{1}.start+...
                    fields{behavior.events.trialConditions(t)}{cellIDs(cell)}{1}.stop)./2;
            else
                com(cell) = nan;
                center(cell) =nan;
            end
        end
%         com(abs(com-1)<10)=nan;
%         com(abs(center-length(behavior.events.map{behavior.events.trialConditions(t)}.x))<10)=nan;
%         center(abs(center-1)<10)=nan;
%         center(abs(center-length(behavior.events.map{behavior.events.trialConditions(t)}.x))<10)=nan;
        predicted_pos(tt) = nanmean(center);
        clear com center;
        else
            predicted_pos(tt) = nan;
        end
           
%        subplot(3,2,2)
%        plot(spikes.spindices((InIntervals(spikes.spindices(:,1),thetaCycleQuartilesTS(trial_thetaCycleQuartilesTS(tt),[1 5]))),1)-thetaCycleQuartilesTS(trial_thetaCycleQuartilesTS(tt),[3]),predicted_lastquart(tt)-actual_pos(tt),'.k')
%        hold on   
%        ylabel('predict last quart - actual')
%        xlabel('position of spike relative to theta trough')
       
   end 
   subplot(3,2,1)
   plot(actual_vel,predicted_lastquart-actual_pos,'.k')
   hold on
   title('actual vel vs  actual - pred lastquart ("lookahead")')
   xlabel('actual vel')
   ylabel('pred lastquart - actual')

   subplot(3,2,3)
   plot(actual_pos, predicted_lastquart-actual_pos,'.k')
   hold on;
   title('actual vs actual - pred lastquart ("lookahead")')
   xlabel('actual')
   ylabel('pred lastquart - actual')
   
   subplot(3,2,4)
   plot(actual_pos, predicted_pos-actual_pos,'.k')
   hold on;
   title('actual vs predicted thetaCycleQuartilesBins 1-3')
   xlabel('actual')
   ylabel('pred thetaCycleQuartilesBins 1-3 - actual')
   
   
   subplot(3,2,5)
   plot(actual_pos,predicted_pos,'.k')
   hold on
   title('pos vs pred pos first thetaCycleQuartilesBins')
   xlabel('actual')
   ylabel('predicted')
   
   subplot(3,2,6);
   plot(actual_pos,predicted_lastquart,'.k')
   hold on
   title('pos  vs pred lastquart')
   xlabel('actual')
   ylabel('predicted')
   
   pause(.001)
   clear pred* actual*
   end
end


end





