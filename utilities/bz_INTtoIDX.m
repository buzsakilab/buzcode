function [ IDX,timestamps ] = bz_INTtoIDX(INT,varargin)
%[IDX] = bz_INTtoIDX(INT) Converts state on/offsets to vector of indices
%
%INPUT
%   INT:    {nstates} cell array of [nintervals x 2] start and end times.
%                       (optional) can be TSObject intervalSet
%           -or-
%           structure with fields: INT.NAMEstate
%   (optional)
%   'statenames'    cell array of state names (position correponds to state number)
%   'sf'            desired sampling frequency of the output vector (default: 1s)
%   'length'        desired length of the index vector (default: max end time)
%
%OUTPUT
%   IDX:    [len x 1] vector of state indices, where states are identified by
%           integers starting from 1, 0 are unmarked.
%           -or-
%           structure with fields:  IDX.statenames
%           (see buzcode wiki)      IDX.states
%                                   IDX.timestamps
%   timestamps
%
%DLevenstein 2015
%Updated 2018 for buzcode
%% DEV
%INT = SleepState.ints;
%statenames = {'WAKE','','NREM','','REM'};
%%
p = inputParser;
addParameter(p,'statenames',[],@iscell)
addParameter(p,'sf',1)
addParameter(p,'length',Inf)
parse(p,varargin{:})
statenames = p.Results.statenames; 
sf = p.Results.sf; 
len = p.Results.length; 
%% Deal with Input Types
%For Buzcode ints structure
if isstruct(INT)
    STRUCTIN = true;
    %Get the NAMEs of the states from the structure (INT.NAMEstate)
    fields = fieldnames(INT);
    endState = cellfun(@(X) contains(X,'state'),fields)';
    fieldstates = cellfun(@(X) char(extractBefore(X,'state')),fields(endState),'uniformoutput',false)';

    if ~isempty(statenames) %Check if there are any states the user didn't put in
        samestates = ismember(fieldstates,statenames);
        %Here: if you find any states in struct not in input, prompt user
        for ss = 1:length(samestates); if samestates(ss)==0
                statenum = inputdlg(['What number would you like state ''',fieldstates{ss},''' to be']);
                statenum = str2num(statenum{1});
                statenames{statenum} = fieldstates{ss};
        end; end     
    else
        statenames = fieldstates;
    end
    
    %Convert to the cell array format needed for this function 
    for ss = 1:length(statenames)
        if isempty(statenames{ss})
            continue; 
        elseif  ~ismember(statenames{ss},fieldstates) 
            INTtemp{ss} = [];
            continue; 
        end
        INTtemp{ss} = INT.([statenames{ss},'state']);
    end
    INT = INTtemp;
else
    STRUCTIN = false;
end

%TSToolbox
if isa(INT,'intervalSet')
    INT = {[Start(INT,'s'), End(INT,'s')]};
end

%%

%Convert from seconds to dt = 1/sf
INT = cellfun(@(X) X*sf,INT,'UniformOutput',false);

%Getting the length of the output
if isinf(len)
    allints = cat(1,INT{:});
    len = ceil(max(allints(:,2)));
end

%%

IDX = zeros(len,1);

numstates = length(INT);
statenums = 1:numstates; %for possible later implementation of non 1:numstates nums
for ss = statenums
    if isempty(INT{ss}); continue; end
    stateints = round(INT{ss});
    
    stateints(stateints==0)=1;
    stateints(isinf(stateints))=len;
    
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

timestamps = [1:length(IDX)]'./sf;

%% Convert to structure
if STRUCTIN
    IDXstruct.states = IDX;
    IDXstruct.timestamps = timestamps;
    IDXstruct.statenames = statenames;
    IDX = IDXstruct;
end

end
