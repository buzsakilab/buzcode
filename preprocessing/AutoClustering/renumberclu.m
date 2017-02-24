function clu=renumberclu(clu)

ix = unique(clu);
for ii=1:length(ix)
    clu(clu==ix(ii))=ii;
end
