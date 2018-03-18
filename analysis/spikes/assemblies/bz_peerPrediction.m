function [dev devControl] =...
    bz_peerPrediction(spikeTimes,winRange,extraPredictors,pairsToRun)
% INPUT
% 
%     spikeTimes - (M x N x D) matrix of spike times. M is the cell number, 
%                  N is the trial number, and D is length of a trial in
%                  milliseconds
%     winRange - vector that represents the range of temporal smoothing 
%                 windows over which to run the assembly analysis 
%                 (default is 0:500 ms)
%     extraPredictors - M x D matrix of variables to regress out (i.e.
%                       position, velocity, theta phase)
%     pairsToRun - (P x 2) matrix of cell pairs that will be analyzed 
%                  (default analyzes all cell pairs within spikeTimes matrix)
%
% OUTPUT
%      dev - (W x M x N) matrix of deviance from solution vectors for each 
%              cell for each trial at each temporally smoothed window size
%      stats - cell array of statistics returns by glmfit.m at each 
%                temporally smoothed window size
%      *Control - same variables as above with trial order randomly shuffled
%
% analysis used in Harris 2003 and Tingley 2015 for finding and quantifying
% cell assemblies.  This code uses a GLM to predict the spike train of a
% cell given a peer's spike train.  
%
% david tingley 2015




if nargin < 1
    error('Must have at least one argument: cell_assembly(spikeTimes)')
elseif nargin == 1
    numCells = size(spikeTimes,1);
    winRange = 0:500;
    pairsToRun = [];
    for n = 1:numCells
        for nn = n:numCells
            if n ~= nn
        pairsToRun = [pairsToRun; n,nn];
            end
        end
    end
elseif nargin == 2
    numCells = size(spikeTimes,1);
    pairsToRun = [];
    for n = 1:numCells
        for nn = n:numCells
            if n ~= nn
        pairsToRun = [pairsToRun; n,nn];
            end
        end
    end
end


numTrials = size(spikeTimes,2);
if length(size(spikeTimes))==2
   numTrials = 1;
   spikes(:,1,:) = spikeTimes; 
   spikeTimes = spikes;
end
% EDIT THIS to take into acount matrices with a single trial (M x D) rather
% than (M x N x D)

% NOTE: Indeces in all outputs match the intervals of winRange, which is not
%       necessarily intervals of 1.
c = 1;

for win = winRange
%     tic
    %              stats{c}   devControl(c,:,:) statsControl{c}
       [dev(win+1,:,:) devControl(win+1,:,:)] = ...
           parGLMRun(win,spikeTimes,extraPredictors,pairsToRun,numTrials);
%        c = c+1
       win
%        toc
end
end

%               stats devControl statsControl
function [dev devControl] = ...
    parGLMRun(win,spikeTimes,extraPredictors,pairsToRun,numTrials)
warning off;
pred_last = 0;
for pair = 1:size(pairsToRun(:,1),1)
%     tic
pred = pairsToRun(pair,1);
act = pairsToRun(pair,2);
if sum(spikeTimes(act,:,:))+sum(spikeTimes(pred,:,:)) > 100
   %% create temporally smoothed predictor spike trains here
if pred_last ~= pred
    for temp = 1:size(spikeTimes,3)
        if temp>win & temp+win<size(spikeTimes,3)
            if win == 0 
            smoothedTrains(:,temp)= spikeTimes(pred,:,temp);
            elseif win > 0
            smoothedTrains(:,temp)=  sum(spikeTimes(pred,:,(temp-win:temp+win))/(win*2),3);
            end
        end
        if temp<=win
            j = temp-1;
            smoothedTrains(:,temp)=  sum(spikeTimes(pred,:,(temp-j:temp+win))/(length(1:temp+win)),3);
        end
        if temp+win>=size(spikeTimes,3)
            j = size(spikeTimes,3)-temp;
            smoothedTrains(:,temp)=  sum(spikeTimes(pred,:,(temp-win:temp+j))/(length(temp-win:size(spikeTimes,3))),3);
        end
%         smoothedTrains(:,temp) = smooth(spikeTimes(pred,:,temp),win);
    end
end
         %% Run GLM on individual trials/instances
        randOrder = randperm(numTrials); 
        for trial = 1:numTrials

            actual = double(squeeze(spikeTimes(act,trial,:))); % GLM's don't like singles...?
            predictor = double(smoothedTrains(trial,:));
            predictorControl = double(smoothedTrains(randOrder(trial),:));
            
%                                stats(pair,trial)
            [results dev(pair,trial) ] = ...
            glmfit([predictor;extraPredictors]',actual,'normal');
%             yhat(pair,trial,:) = glmval(results,[predictor; 1:size(spikeTimes,3)]','identity');

            if numTrials == 1 % we need another way to shuffle if only one trial is given
                for iter = 1:5
                    predictorControlShifted = circshift(predictorControl,...
                        round(rand*length(predictorControl)),2);
                    [resultsControl(:,iter) devControl(pair,trial,iter)] = ...
                    glmfit([predictorControlShifted;extraPredictors]',actual,'normal');
                end
            else
                [resultsControl devControl(pair,trial)] = ...
                    glmfit([predictorControl;extraPredictors]',actual,'normal');
            end
            
%             yhatControl(pair,trial,:) = glmval(results,[predictorControl; 1:size(spikeTimes,3)]','identity');            
        end   
     pred_last = pred;
%      toc
 end
end
end
