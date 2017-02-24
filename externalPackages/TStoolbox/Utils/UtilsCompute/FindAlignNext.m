function ix = FindAlignNext(D, t)



nIX = length(t);

ix = zeros(size(t));
for iIX = 1:nIX
   ix(iIX) = AlignNext(D, t(iIX));
end
