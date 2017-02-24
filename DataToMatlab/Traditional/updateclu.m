function clu = updateclu(clu,ix,newclu)

for ii=1:length(ix)
    clu(clu==ix(ii))=newclu;
end

% ix = unique(clu);
% for ii=1:length(ix)
%     clu(clu==ix(ii))=ii;
% end
