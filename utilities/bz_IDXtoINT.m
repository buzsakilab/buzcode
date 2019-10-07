function [ INT ] = bz_IDXtoINT( IDX ,varargin)
%bz_IDXtoINT(IDX) Converts state indices to state on/offsets
%
%INPUT
%   IDX:    structure with fields:  IDX.statenames
%           (see buzcode wiki)      IDX.states
%                                   IDX.timestamps
%           -or-
%           [len x 1] vector of state indices, where states are identified by
%           integers starting from 1, 0 are unmarked. 
%           If IDX is intered as a vector, 'timestamps' and 'statenames' 
%           can be entered as optional inputs. Alternatively,
%           Timestamps will be assumed to be [1:length(IDX)]*dt and states
%           will be named for you
%   
%   (options)
%   'numstates' (optional)  number of interval types (for use
%   'statenames'
%   'timestamps'
%   'nameStates'
%   'dt'
%   'jumptol'   tolernce for jumps, in units of dt (default: 2). If
%               adjacent timestamps are bigger than the tolerence, will
%               end/start interval around the jump
%   
%OUTPUT
%   INT:    {nstates} cell array of intervals - start and end times
%           -or-
%           structure with fields: INT.NAMEstate
%
%DLevenstein 2015-16
%%
p = inputParser;
addParameter(p,'statenames',[],@iscell)
addParameter(p,'timestamps',[])
addParameter(p,'numstates',[])
addParameter(p,'dt',[])
addParameter(p,'nameStates',false)
addParameter(p,'jumptol',2)
parse(p,varargin{:})
statenames = p.Results.statenames; 
dt = p.Results.dt; 
numstates = p.Results.numstates; 
nameStates = p.Results.nameStates; 
jumptol = p.Results.jumptol; 

%%
if isstruct(IDX)
    if isfield(IDX,'statenames')
        statenames = IDX.statenames;
    end
    timestamps = IDX.timestamps;
    IDX = IDX.states;
end


if islogical(IDX)
    IDX = double(IDX); 
end

if exist('statenames','var')
    numstates = length(statenames);
elseif isempty('numstates')
    numstates = max(IDX);
end

if isempty(statenames)
    statenames = cell(1,numstates);
end


states = 1:numstates;

if isrow(IDX)
    IDX = IDX';
end
%%

if isempty(dt)
   dt = mode(diff(timestamps));
end
%For timestamps with breaks Fill in with state 0 and timestamp nan 
if any(diff(timestamps)>(jumptol.*dt))
    jumps = find(diff(timestamps)>(jumptol.*dt));
    for jj = length(jumps):-1:1
        IDX = [IDX(1:jumps(jj));0;IDX(jumps(jj)+1:end)];
        timestamps = [timestamps(1:jumps(jj));nan;timestamps(jumps(jj)+1:end)];
    end
   %
end

IDX = [0; IDX; 0]; %So that diff will put the index on the right timestamp
%%
for ss = 1:numstates
    statetimes = IDX==states(ss);
    %Get indices of state on/offsets
    stateints = [find(diff(statetimes)==1) find(diff(statetimes)==-1)-1];
    stateints = [timestamps(stateints(:,1)) timestamps(stateints(:,2))]; %Get timestamps of state on/offsets
    
    if isempty(stateints) && isempty(statenames{ss});continue;end
    
    %Get the name
    if ~isempty(statenames{ss})
        newname = strcat(statenames{ss},'state');
    else
        if nameStates
            usernamedstate = inputdlg(['What is the name of state ',num2str(ss),'?']);
            statenames{ss} = usernamedstate{1};
            newname = strcat(statenames{ss},'state');
        else
            newname = strcat('state',ss);
        end
    end
    
    INT.(newname) = stateints;
    clear stateints
end



end
