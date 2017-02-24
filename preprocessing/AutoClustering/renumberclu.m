function [clu log] =renumberclu(clu,log)
% re-orders clu to start from index 1 and iterate through the number of
% clusters

% rewritten by David Tingley, 2017

cluster_names = unique(clu);
for i=1:length(cluster_names)
    if cluster_names(i) ~= 0  % keep cluster 0 (noise) as 0
    clu(find(clu==cluster_names(i))) = i;
     log = [log sprintf('%d -> %d; reordering clusters\n',cluster_names(i),i)];
    end
end

