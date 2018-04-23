function [ IDX ] = INTtoIDX_ss(INT,len,sf)
%[IDX] = INTtoIDX(INT,len,sf) Converts state on/offsets to vector of indices
%
%INPUT
%   INT:    {nstates} cell array of [nintervals x 2] start and end times.
%                       (optional) can be TSObject intervalSet
%   len:    length of index vector
%   sf:     desired sampling frequency of the output vector
%
%OUTPUT
%   IDX:    [len x 1] vector of state indices, where states are identified by
%           integers starting from 1, 0 are unmarked.
%
%Last Updated: 11/15/15
%DLevenstein
%%
if isa(INT,'intervalSet')
    INT = {[Start(INT,'s'), End(INT,'s')]};
end

if exist('sf','var')
    INT = cellfun(@(X) X*sf,INT,'UniformOutput',false);
end

if ~exist('len','var')
    allints = cat(1,INT{:});
    len = max(allints(:));
end

IDX = zeros(len,1);

numstates = length(INT);
for ss = 1:numstates
    stateints = INT{ss};
    numints = length(stateints(:,1));
    for ii = 1:numints
        IDX(stateints(ii,1):stateints(ii,2))=ss;
    end
end

switch numstates
    case 1
        IDX = logical(IDX);
    otherwise
end

if length(IDX)>len
    display("your states extend past your len...")
    IDX = IDX(1:len);
end


end

