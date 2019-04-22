function [] = bz_compareReplay(spikes,ripples,template,include)


% spikes  - buzcode cellinfo file (requires spikes.times, an Nx1 cell
%           array of timestamps in seconds for each neuron; spikes.UID and spikes.spindices)
% ripples - buzcode ripples events file (only requires ripples.timestamps,
%           an Nx2 matrix of start/stop times for each event)
% template - NxD matrix of N cells and D positions, average firing rates
% include - indices (1:N) of cells (place cells) to keep

warning off
binSize = .01;
overlap = 1;
spkmat = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize * overlap);


if ~isempty(include)
    keep = intersect(include,find(sum(template')>0)); % uncomment to keep all HPC cells that fired
else
    keep = find(sum(template')>0); % just remove cells w/ zeros spikes in template
end

% normalize template
% for i = 1:size(template,1)
%     template(i,:) = mean_norm(template(i,:));
% end
            

for event = 1:size(ripples.timestamps,1)
                % start bigger than the event itself
%                 start = round((round(ripples.timestamps(event,1) * 1000)) ./ (spkmat.dt*1000));
%                 stop = round((round(ripples.timestamps(event,2) * 1000)) ./ (spkmat.dt*1000));
                start = round((round(ripples.peaks(event) * 1000)-50) ./ (spkmat.dt*1000));
                stop = round((round(ripples.peaks(event) * 1000)+50) ./ (spkmat.dt*1000));
                
               
                    for spk = 1:size(spkmat.data,2)
                            data(:,spk) = (spkmat.data(start:stop,spk)')';  
                            counts(:,spk) = (spkmat.data(start:stop,spk)')';  
                    end

                    % cut 0 rate bins from the start/end
                    while sum(counts(1,keep)) < 0 & size(counts,1) > 1
                       data = data(2:end,:);
                       counts = counts(2:end,:);
                       start = start + 1;
                    end
                    while sum(counts(end,keep)) < 0 & size(counts,1) > 1
                       data = data(1:end-1,:);
                       counts = counts(1:end-1,:);
                       stop = stop-1;
                    end
                        

                if stop < size(spkmat.data,1) & size(data,1) > 4 & sum(sum(counts(:,keep))>0) > 4
                    % norm data matrix
                    data = data ./ binSize;
                    % calc the posterior matrix
                    [Pr, prMax] = placeBayes(data(:,keep), template(keep,:), spkmat.dt);
                    
                    % now calc correlations w/ the data...
                    [bayesRankOrder(event) ]= corr([1:length(prMax)]',prMax,'rows','complete');
                    % linear weighted correlation
                    [linearWeighted(event) outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));
                    
                    % extra info to capture about event
                    nCells(event) = sum(sum(counts(:,keep))>0);
                    nSpks(event) = sum(sum(counts(:,keep)));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rank order calculation
                    % get all spks from event
%                   spks = find(InIntervals(spikes.spindices(:,1),[ripples.timestamps(event,1) ripples.timestamps(event,2)]));
%                     spks = find(InIntervals(spikes.spindices(:,1),[spkmat.timestamps(start) spkmat.timestamps(stop)]));
%                     
%                     % for each cell get first and avg spk times
%                     UIDs = spikes.UID;
%                     for cell=1:size(data,2)
%                         ts((cell)) = nanmean(spikes.spindices(spikes.spindices(spks,2)==UIDs(cell),1));
%                         temp = spikes.spindices(spikes.spindices(spks,2)==UIDs(cell),1);
%                         if ~isempty(temp)
%                             ts_first((cell)) = temp(1);
%                         else
%                             ts_first((cell)) = nan;
%                         end
%                     end
% 
%                     idx = intersect(find(~isnan(ts)),keep); % only take cells that spiked...               
%                     [a b ord_template] = sort_cells(template(idx,:));
%                     [a ord_avg] = sortrows(ts(idx)','descend');
%                     [a ord_firstSpk] = sortrows(ts_first(idx)','descend');
%                     clear ts 
                    idx = intersect(find(sum(data)>0),keep);
                    [a b ord_template] = sort_cells(template(idx,:));
                    [a ord_firstSpk] = sortrows(data(:,idx)','descend');
                    [rankOrder(event) pvals(event)] = corr(ord_template,ord_firstSpk,'rows','complete');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    
                  
                    if sum(sum(spkmat.data(start:stop,:)))> 5 * overlap & sum(~isnan(sum(Pr')))>5
                        [slope_hpc(event) integral_hpc(event)] = Pr2Radon(Pr');
                    else 
                        bayesRankOrder(event) = nan;
                        bayesRankOrder_shuf(event,1:100) = nan;
                        rankOrder(event) = nan;
                        linearWeighted(event) = nan;
                        linearWeighted_shuf(event,1:100) = nan;
                        slope_hpc(event) = nan;
                        integral_hpc(event) =nan;
                        slope_hpc_shuf(event) = nan;
                        integral_hpc_shuf(event,1:100) = nan;
                    end
                else
                    bayesRankOrder(event) = nan;
                    bayesRankOrder_shuf(event,1:100) = nan;
                    rankOrder(event) = nan;
                    linearWeighted(event) = nan;
                    linearWeighted_shuf(event,1:100) = nan;
                    slope_hpc(event) = nan;
                    integral_hpc(event) =nan;
                    slope_hpc_shuf(event) = nan;
                    integral_hpc_shuf(event,1:100) = nan;
                    Pr = [];
                end

                nCells(event) = length(keep);
                spkCount(event) = sum(sum(data(:,keep)))./overlap;
                eventDuration(event) = (stop-start)*spkmat.dt;
           
                


subplot(4,2,1)
scatter(rankOrder,linearWeighted,'.k')
title('rank ord VS linear weighted')

subplot(4,2,2)
scatter(linearWeighted,bayesRankOrder,'.k')
title('linear weighted VS bayesian rank ord')
subplot(4,2,3)
scatter(rankOrder,integral_hpc,'.k')
title('rank ord VS radon in')

subplot(4,2,4)
scatter(linearWeighted,integral_hpc,'.k')
title('linear weighted VS radon int')

subplot(4,2,5)
scatter(rankOrder,slope_hpc,'.k')
title('rank ord VS radon slope')

subplot(4,2,6)
scatter(linearWeighted,slope_hpc,'.k')
title('linear weighted VS radon slope')

subplot(4,2,7)
imagesc(data')
title(rankOrder(end))

subplot(4,2,8)
imagesc(Pr')
title(linearWeighted(end))



pause(.01)
clear data counts;
event


end


data.linearWeighted = linearWeighted;
