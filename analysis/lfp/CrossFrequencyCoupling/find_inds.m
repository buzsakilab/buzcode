function inds = find_inds(targets, ref_vec, PAR)
% inds = find_inds(targets, ref_vec)
%
% Created by: Erik Schomburg, 2011

if (nargin < 3)
    PAR = 0;
end

inds = zeros(size(targets));
if PAR
    parfor i=1:length(targets(:))
        [~,inds(i)] = min(abs(targets(i) - ref_vec));
    end
else
    for i=1:length(targets(:))
        [~,inds(i)] = min(abs(targets(i) - ref_vec));
    end
end
