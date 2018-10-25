function [ spiketimes_jitt ] = JitterSpiketimes(spiketimes,jitterwin,varargin)
%[ spiketimes ] = JitterSpiketimes(spiketimes,jitterwin) jitters spiketimes
%of all cells within a designated window.
%
%INPUT
%   spiketimes  {Ncells} cell array of [Nspikes] vectors of spiketimes for
%               each cell (seconds)
%   jitterwin   time window within which to jitter spike times (seconds)
%   
%   (options)
%       'jittertype'    'window' or 'centered'. 
%                       -'centered' jitters each spike around its time.
%                       -'window' jitters within its window (more rigorous, default)
%                           see Jonathan Platkewitz paper
%
%OUTPUT
%   spiketimes  jittered spiketimes
%
%
%DLevenstein 2016
%%
%parse args
p = inputParser;
addParameter(p,'jittertype','window')

parse(p,varargin{:})
jittertype = p.Results.jittertype;


%%
numcells = length(spiketimes);
if isa(spiketimes,'tsdArray')
    for c = 1:numcells
        spiketimestemp{c} = Range(spiketimes{c},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
end


%Do this with cellfun....
switch jittertype
    case 'centered'
        for c = 1:numcells
            spiketimes_jitt{c} = spiketimes{c}+2*jitterwin*rand(size(spiketimes{c}))-jitterwin;
        end
    case 'window'
        %Question: add global offset, so not always starting windows at same
        %points in recording?
        spiketimes_jitt = cellfun(@(X) ...
            jitterwin.*(floor(X/jitterwin)+rand(size(X))),...
            spiketimes,'UniformOutput',false);       
end

%Re-sort
spiketimes_jitt = cellfun(@(X) sort(X),spiketimes_jitt,'UniformOutput',false);

end

