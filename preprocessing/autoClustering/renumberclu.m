function [new_clu log] =renumberclu(clu,log)
% re-orders clu to start from index 1 and iterate through the number of
% clusters

% rewritten by David Tingley, 2017
new_clu = zeros(length(clu),1);

cluster_names = unique(clu);
for i=1:length(cluster_names)
    if cluster_names(i) ~= 0  % keep cluster 0 (noise) as 0
    new_clu(find(clu==cluster_names(i))) = i+2;
     log = [log sprintf('%d -> %d; reordering clusters\n',cluster_names(i),i+2)];
    end
end

