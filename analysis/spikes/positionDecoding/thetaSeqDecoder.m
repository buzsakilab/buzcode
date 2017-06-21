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
%   trialIntervals
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

numSpikesThresh = 3;

%% first let's get the ratemaps for these cells, and find the place field 'centers'

disp('getting ratemaps and place field statistics..')

[rateMap countMap occuMap phaseMap] = spaceRateMap(spikes.times,pos,map,mapping,trials,lfp);
% remove high firing rate neurons
for i=1:length(spikes.times)
    FR = length(spikes.times{i}) ./ lfp.duration;
    if FR > 5
        spikes.times{i} = [];
    end
end
[spktimes spkIDs] = spikes2sorted(spikes.times);

for i=1:8
for j=1:length(trials{i})
    if size(trials{i}{j},2) == 5
        t_ints{i}(j,:)= [trials{i}{j}(1,5);trials{i}{j}(end,5); i; j];
    else
        t_ints{i}(j,:)= [trials{i}{j}(1,1);trials{i}{j}(end,1); i; j];
    end
end
end
trialIntervals = cell2mat(t_ints');

% let's calculate the population rate
pop_rate = zeros(length(lfp.data),1);
for i=1:length(spikes.times)
    pop_rate(round(1250*spikes.times{i})) = pop_rate(round(1250*spikes.times{i})) + 1;
end

% find place fields and their COM's
for i=1:length(rateMap)
    fields{i} = bz_getPlaceFields1D(rateMap{i});
end

%% now let's find all theta cycles and split them into quartiles

[b a] = butter(4,[5/(lfp.samplingRate/2) 16/(lfp.samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
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

for t = 1:size(trialIntervals,1)
   trial_cycles = find(InIntervals(cycles(:,3),trialIntervals(t,1:2))); 
   if ~isempty(trial_cycles)
   for tt = 1:length(trial_cycles)
    	cellIDs = unique(spkIDs((InIntervals(spktimes,cycles(trial_cycles(tt),4:5)))));
        if length(cellIDs) > numSpikesThresh
        for cell = 1:length(cellIDs)
            if ~isempty(fields{trialIntervals(t,3)}{cellIDs(cell)}) & length(fields{trialIntervals(t,3)}{cellIDs(cell)}) == 1
                com(cell) = fields{trialIntervals(t,3)}{cellIDs(cell)}{1}.COM;
            else
                com(cell) = nan;
            end
        end
        predicted_lastquart(tt) = nanmedian(com);
        if size(mapping{trialIntervals(t,3)}{trialIntervals(t,4)},2) == 14
            [a b] = min(abs(mapping{trialIntervals(t,3)}{trialIntervals(t,4)}(:,14) - cycles(trial_cycles(tt),3)));
            actual_pos(tt) = mapping{trialIntervals(t,3)}{trialIntervals(t,4)}(b,13);
        else
            [a b] = min(abs(mapping{trialIntervals(t,3)}{trialIntervals(t,4)}(:,6) - cycles(trial_cycles(tt),3)));
             actual_pos(tt) = mapping{trialIntervals(t,3)}{trialIntervals(t,4)}(b,5);
        end
        
        clear com;
        else
            actual_pos(tt)=nan;
            predicted_lastquart(tt) = nan;
        end
        
        cellIDs = unique(spkIDs((InIntervals(spktimes,cycles(trial_cycles(tt),[1 4])))));
        if length(cellIDs) > numSpikesThresh
        for cell = 1:length(cellIDs)
            if ~isempty(fields{trialIntervals(t,3)}{cellIDs(cell)}) & length(fields{trialIntervals(t,3)}{cellIDs(cell)}) == 1
                com(cell) = fields{trialIntervals(t,3)}{cellIDs(cell)}{1}.COM;
            else
                com(cell) = nan;
            end
        end
        predicted_pos(tt) = nanmedian(com);
        clear com;
        else
            predicted_pos(tt) = nan;
        end
   end 
   subplot(3,2,1)
   plot(actual_pos,predicted_pos-predicted_lastquart,'.k')
   hold on
   title('actual pos vs predicted diff')
   xlabel('actual')
   ylabel('predicted firstquarts - lastquart')
   
   subplot(3,2,3)
   plot(actual_pos, actual_pos-predicted_lastquart,'.k')
   hold on;
   title('actual vs actual - pred lastquart ("lookahead")')
   xlabel('actual')
   ylabel('actual - pred lastquart')
   
   subplot(3,2,4)
   plot(actual_pos, actual_pos-predicted_pos,'.k')
   hold on;
   title('actual vs predicted quarts 1-3 ')
   xlabel('actual')
   ylabel('actual - pred quarts 1-3')
   
   
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





