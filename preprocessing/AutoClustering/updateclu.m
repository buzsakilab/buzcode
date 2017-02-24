function clu = updateclu(clu,ix,newclu)
% updates identities in parts of a clu file

%% INPUTS

% clu - N x 1 vector of cluster identities for each extracted waveform
% ix - M x 1 vector of indices where new identities will be assigned
% newclu - integer value for the new cluster identity

%% OUTPUTS
% clu - N x 1 vector of updated cluster identities

% updated 2017/02/24 by David Tingley

for ii=1:length(ix)
    clu(clu==ix(ii))=newclu;
end

% ix = unique(clu);
% for ii=1:length(ix)
%     clu(clu==ix(ii))=ii;
% end
