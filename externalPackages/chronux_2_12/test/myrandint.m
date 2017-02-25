function ranInt = myrandint(outputRow,outputCol,outputRange,varargin)
% MYRANDINT(M,N,RANGE) is an M-by-N matrix with random integer entries 
% drawn with replacement from elements of vector RANGE.  The elements in
% vector RANGE do not need to be contiguous or unique.  (Actually, they do
% not even need to be integers: The function works the exact same way with
% noninteger elements, but a warning is generated to alert the user that 
% noninteger elements are being sampled.) 
%
% To specify a contiguous integer range from Xlow to Xhi, use RANGE = [Xlow:Xhi].  
% 
% MYRANDINT(M,N,RANGE,'noreplace') is an M-by-N matrix with random integers 
% drawn without replacement.
%
% This function is based around RAND and RANDPERM, and is intended as a
% modest imitation of Comm Toolbox's RANDINT.


if isequal(size(outputRange),[1 2]) && ~isequal(outputRange(1),outputRange(2)-1),
    warning('To specify a range [low high] use [low:high].')
end
if ~isequal(round(outputRange),outputRange),
    warning('Specified RANGE contains noninteger values.')
end
if ~isequal(length(outputRange),length(outputRange(:))),
    error('Range must be a vector of integer values.')
end

numElements = outputRow*outputCol;

if isempty(varargin),
    
    ranInt = zeros(outputRow,outputCol);
    randIx = floor((length(outputRange))*rand(size(ranInt))) + 1;
    ranInt = outputRange(randIx);
    if ~isequal(size(randIx),size(ranInt)),
        ranInt = reshape(ranInt,size(randIx));
    end
    
elseif isequal(varargin{1},'noreplace'),
    
    if numElements > length(outputRange),
        error('Not enough elements in range to sample without replacement.')
    else
        % Generate full range of integers
        XfullShuffle = outputRange(randperm(length(outputRange)));
        % Select the first bunch:
        ranInt = reshape(XfullShuffle(1:numElements),outputRow,outputCol);     
    end    
    
else
    error('Valid argument is ''noreplace''.')
end


