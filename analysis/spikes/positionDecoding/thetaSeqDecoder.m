function [] = thetaSeqDecoder(spikes,lfp,behav)
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
% decoding done in the following paper: Wikenheiser and Redish. Hippocampal
% theta sequences reflect current goals 2015
% 
% written by david tingley, 2017

numSpikesThresh = 2;

%% first let's get the ratemaps for these cells, and find the place field 'centers'

disp('getting ratemaps and place field statistics..')
% if size(pos,2) > 6
% [rateMap countMap occuMap phaseMap] = spaceRateMap(spikes.times,pos,map,mapping,trials,lfp);
% else
% [rateMap countMap occuMap phaseMap] = spaceRateMap_old(spikes.times,pos,map,mapping,trials,[lfp.timestamps double(lfp.data)]);
% end
% spikes = bz_GetSpikes('region','hpc');
for j=1:length(spikes.times)
    for i=1:length(behavior.events.trialConditions)
        if strcmp(behavior.trackingType,'led')
             [curve stats] = FiringCurve([behavior.events.trials{i}.timestamps behavior.events.trials{i}.mapping],spikes.times{j},'nBins',80);
        elseif strcmp(behavior.trackingType,'optitrack')
             [curve stats] = FiringCurve([behavior.events.trials{i}.timestamps behavior.events.trials{i}.mapping],spikes.times{j},'nBins',200);
        end
        rateMap(j,i,:) = curve.rate;
    end
end
% remove high firing rate neurons
for i=1:length(spikes.times)
    FR = length(spikes.times{i}) ./ lfp.duration;
    if FR > 5
        spikes.times{i} = [];
    end
end
[spktimes spkIDs] = spikes2sorted(spikes.times);


% let's calculate the population rate
pop_rate = zeros(length(lfp.data),1);
for i=1:length(spikes.times)
    pop_rate(round(1250*spikes.times{i})) = pop_rate(round(1250*spikes.times{i})) + 1;
end

% find place fields and their COM's
for i=1:length(unique(behavior.events.trialConditions))
    fields{i} = bz_getPlaceFields1D(rateMap(:,behavior.events.trialConditions==i,:),'minPeakRate',3,'maxFieldWidth',25,'minFieldWidth',4);
end

%% now let's find all theta cycles and split them into quartiles
[b a] = butter(4,[1/(lfp.samplingRate/2) 4/(lfp.samplingRate/2)],'bandpass');
delta= FiltFiltM(b,a,double(lfp.data(:,1)));
delta_pow = fastrms(delta,1000);

[b a] = butter(4,[6/(lfp.samplingRate/2) 10/(lfp.samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
pow = fastrms(filt,250);

% phases = angle(hilbert(FiltFiltM(b,a,double(lfp.data(:,1)))));
[blah peaks] = findpeaks(filt);
[blah troughs] = findpeaks(-filt);
% [blah peaks] = findpeaks(cos(phases));
% [blah troughs] = findpeaks(-cos(phases));

disp('segmenting theta cycles into quartiles by peaks/troughs..')
if peaks(1) < troughs(1) 
        peaks(1) = [];
end
    
for i=1:length(peaks) - 1
    quarts(i,:) = [peaks(i) (peaks(i)+troughs(i+1))/2 troughs(i+1)  (peaks(i+1)+troughs(i+1))/2 peaks(i+1)];
end

% check that there are enough spikes for last theta cycle quartile
keep=[];exclude=[];

for i = 1:size(quarts,1)
   if sum(pop_rate(quarts(i,4):quarts(i,5))) > numSpikesThresh
       keep = [keep; i];
   else
       exclude = [exclude;i];
   end
end

quarts = quarts(keep,:);
disp(['searching through ' num2str(size(quarts,1)) ' candidate theta cycles'])

%% now we calculate the distance between the animals position and the COM 
%% for cells spiking in the last quartile of each theta sequence
cycles = quarts ./  lfp.samplingRate;

for t = 1:size(behavior.events.trialIntervals,1)
   trial_cycles = find(InIntervals(cycles(:,3),behavior.events.trialIntervals(t,1:2))); 
   if ~isempty(trial_cycles)
   for tt = 1:length(trial_cycles)
    	cellIDs = (spkIDs((InIntervals(spktimes,cycles(trial_cycles(tt),[4 5])))));
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
        
        [a b] = min(abs(behavior.events.trials{t}.timestamps - cycles(trial_cycles(tt),3)));
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
        
        cellIDs = (spkIDs((InIntervals(spktimes,cycles(trial_cycles(tt),[1 5])))));
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
           
       subplot(3,2,2)
       plot(spktimes((InIntervals(spktimes,cycles(trial_cycles(tt),[1 5]))))-cycles(trial_cycles(tt),[3]),predicted_lastquart(tt)-actual_pos(tt),'.k')
       hold on   
       ylabel('predict last quart - actual')
       xlabel('position of spike relative to theta trough')
       
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
   title('actual vs predicted quarts 1-5 ')
   xlabel('actual')
   ylabel('pred quarts 1-5 - actual')
   
   
   subplot(3,2,5)
   plot(actual_pos,predicted_pos,'.k')
   hold on
   title('pos vs pred pos first quarts')
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





