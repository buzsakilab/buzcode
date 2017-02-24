function [clu log] =renumberclu(clu,log)
% rewritten by David Tingley, 2017

cluster_names = unique(clu);
for i=1:length(cluster_names)
    if cluster_names(i) ~= 0
    clu(find(clu==cluster_names(i))) = i;
     log = [log sprintf('%d -> %d; reordering clusters\n',cluster_names(i),i)];
    end
end

