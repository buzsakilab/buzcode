function [] = bz_compareReplay(spikes,ripples,template,include)


% spikes  - buzcode cellinfo file
% ripples - buzcode ripples events file
% template - NxD matrix of N cells and D positions, average firing rates
% include - indices (1:N) of cells (place cells) to keep

binSize = .01;
overlap = 1;
spkmat = bz_SpktToSpkmat(spikes.times,'overlap',overlap,'binSize',binSize * overlap);

% keep = intersect(find(pf),idx_hpc)-length(idx_ls); % are you a place cell, in HPC?
if ~isempty(include)
    keep = intersect(include,find(sum(template')>0)); % uncomment to keep all HPC cells that fired
else
    keep = find(sum(template')>0); 
end
% keep(sum(rates')==0) = [];


% normalize template
for i = 1:size(template,1)
    template(i,:) = mean_norm(template(i,:)')';
end
            

for event = 1:length(ripples.peaks)
                % start bigger than the event itself
                start = round((round(ripples.timestamps(event,1) * 1000)-50) ./ (spkmat.dt*1000));
                stop = round((round(ripples.timestamps(event,2) * 1000)+50) ./ (spkmat.dt*1000));

                if stop < size(spkmat.data,1) & stop-start > 4
                    for spk = 1:size(spkmat.data,2)
                            data(:,spk) = mean_norm(spkmat.data(start:stop,spk)')';  
                            counts(:,spk) = (spkmat.data(start:stop,spk)')';  
                    end

                    % cut 0 rate bins from the start/end
                    while sum(counts(1,keep)) < 0 & size(counts,1) > 1
                       data = data(2:end,:);
                       counts = counts(2:end,:);
                    end
                    while sum(counts(end,keep)) < 0 & size(counts,1) > 1
                       data = data(1:end-1,:);
                       counts = counts(1:end-1,:);
                    end
                    
                    %% get max rate bin
%                     [a b] = max(sum(data,2));
%     %                 or search from middle...
%     %                 b = round(size(data,1)/2);
%                      % find clipping point at beginning
%                      sta = 0;
%                     while sum(counts(b-sta,keep),2) > 0 & b-sta > 1 % max length is +/-50
%                        sta = sta + 1;
%                     end
%     %                  % find clipping point at beginning
%                      sto = 0;
%                     while sum(counts(b+sto,keep),2) > 0 & sto +b < size(data,1)-1  % max length 500 ms
%                        sto = sto + 1;
%                     end
% 
%                     % redefine start/stop by pop burst..
%                     data = data(b-sta:b+sto,:);
    
    

                if stop < size(spkmat.data,1) & size(data,1) > 4 & sum(sum(counts(:,keep))>0) > 4
                    
                    [Pr, prMax] = placeBayes(data(:,keep), template(keep,:), spkmat.dt);
                    
                    % now calc correlations w/ the data...
                    [bayesRankOrder(event) ]= corr([1:length(prMax)]',prMax,'rows','complete');
                    
                    [linearWeighted(event) outID] = makeBayesWeightedCorr1(Pr,ones(size(Pr,1),1));
                    
                    nCells(event) = sum(sum(counts(:,keep))>0);
                    nSpks(event) = sum(sum(counts(:,keep)));

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    idx = intersect(find(nanmean(data)~=0),keep); % only take cells that spiked...               
                    [a b ord_template] = sort_cells(template(idx,:));
                    [a b ord_firstSpk] = sort_cells(data(:,idx)');
                    for i =1:size(data,2)
                        [ts(i)]=mean(find(data(:,i)>0));
                    end
                    [a b ord_avg] = sort_cells(ts(idx)'); clear ts
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
                end
                pause(.01)
                clear data counts;
                event


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


end

end
