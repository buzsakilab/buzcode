function R = SelectAlongFirstDimension(IN, f)

% R = SelectAlongFirstDimension(IN,f)
%
% equivalent of IN(f) for 1D, IN(f,:) for 2D, IN(f,:,:) for 3D, etc...
%
% INPUTS:
%    IN = any matrix input
%    f = selection indices (such as returned by find)
%
% OUTPUTS:
%    R = same type matrix as IN 
%

% ADR 1998
% version L5.0
% status PROMOTED
% v 5.0 now can handle NaN indices

% new version Francesco P. Battaglia 2004
% now can handle empty matrices 
  
if isempty(IN)
  R = [];
  return
end

  
sz = size(IN);
dim1 = sz(1);
dimRest = sz(2:length(sz));

tmpMatrix = double(reshape(IN, [dim1 prod(dimRest)]));

nonnan = find(~isnan(f));
yesnan = find(isnan(f));
f1 = f;
f1(yesnan) = 1;
tmpMatrix = tmpMatrix(f1,:);
tmpMatrix(yesnan,:) = NaN;

R = reshape(tmpMatrix, [length(f) dimRest]);