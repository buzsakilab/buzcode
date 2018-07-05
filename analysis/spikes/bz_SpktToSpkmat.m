function [spikemat] = bz_SpktToSpkmat(spikes, varargin)
%spikemat = SpktToSpkmat(spiketimes,<options>) takes a 
% 1 x N_neurons cell array of spiketimes  and converts into a t/dt x N spike
% matrix.
%
%INPUTS
%   spikes      a structure with fields spikes.times
%               can also take a {Ncells} cell array of [spiketimes]
%               or [spiketimes, UID] pairs
%
%       (options)
%       'win'       [start stop] time interval to get spike matrix
%       'binsize'   size of your time bins, in seconds
%       'overlap'   overlap of your time bins
%       'dt'  add (use instead of binsize/overlap)
%
%OUTPUT
%
%Return:
%Spike Matrix
%time vector
%Spike time/cell indices for plotting
%
%To Do:
%   -Remove for loop... don't need to go through structure.
%   -T is just silly... make this a reasonable time window able to select
%   only spikes within a given time
%
%
%DLevenstein 2015. Updated 2018 for buzcode
%NOTE: in progres...
%% Options
p = inputParser;
addParameter(p,'win',[]);
addParameter(p,'binsize',0.1);
addParameter(p,'overlap',1);


parse(p,varargin{:})

win = p.Results.win;
binsize = p.Results.binsize;
overlap = p.Results.overlap;

dt = binsize./overlap;



%% Deal With Input Type Variability

%If spiketimes is in a buzcode structure
if isstruct(spikes)
    spiketimes = spikes.times;
else
    spiketimes = spikes;
end

%If spike times is in the form [spiketimes(:,1) cellnum(:,2)], convert to
%cell array
if isa(spiketimes,'numeric') && size(spiketimes,2)==2
    cellnums = unique(spiketimes(:,2));
    for cc = cellnums'
        spiketimestemp{cc} = spiketimes(spiketimes(:,2)==cc,1);
    end
    spiketimes = spiketimestemp;
end

%Take stock of the cells - if there are no cells that's silly, but doesn't
%break.
numcells = length(spiketimes);
if numcells == 0
    spkmat=[];t=[];spindices=[];
    return
end

%Spiketimes can be: tsdArray of cells, cell array of cells, cell array of
%tsdArrays (multiple populations)
if isa(spiketimes,'tsdArray')
    numcells = length(spiketimes);
    for cc = 1:numcells
        spiketimestemp{cc} = Range(spiketimes{cc},'s');
    end
    spiketimes = spiketimestemp;
    clear spiketimestemp
elseif isa(spiketimes,'cell') && isa(spiketimes{1},'tsdArray')
    numpop = length(spiketimes);
    lastpopnum = 0;
    for pp = 1:numpop
        if length(spiketimes{pp})==0
            spiketimes{pp} = {};
            popcellind{pp} = [];
            continue
        end
        for cc = 1:length(spiketimes{pp})
            spiketimestemp{cc} = Range(spiketimes{pp}{cc},'s');
        end
        spiketimes{pp} = spiketimestemp;
        popcellind{pp} = [1:length(spiketimes{pp})]+lastpopnum;
        lastpopnum = popcellind{pp}(end);
        clear spiketimestemp
    end
    spiketimes = cat(2,spiketimes{:});
    numcells = length(spiketimes);
    subpop = 'done';
end

%Time Window
if isempty(win) || isequal(win,[0 Inf])
    t_start = 0; t_end = max(vertcat(spiketimes{:}));
elseif  length(win) == 2
    t_start = win(1); t_end = win(2);
end

%% The Meat of the function

numts = ceil((t_end-t_start)/dt);

%Remove spikes after t_end and before t_start (t_offset+t_start)
spiketimes = cellfun(@(x) x(find(x<t_end)),spiketimes,'UniformOutput',false);
spiketimes = cellfun(@(x) x(find(x>t_start)),spiketimes,'UniformOutput',false);


%Establish Cell Structure... maybe do this with cellfun... or
%cell2struct
cells = cell2struct(spiketimes,'spiketimes',1);
for cell_ind = 1:numcells
    %When Spike? row index for each spike
    cells(cell_ind).spiketimes = cells(cell_ind).spiketimes';
    %Which Cell? column index for each spike
    cells(cell_ind).index4spikes = cell_ind*ones(size(cells(cell_ind).spiketimes));
end

%Make a Spike Matrix
spkmat = zeros(numts,numcells);
%Spike Indices - time
spikes_ind_t = ceil(([cells.spiketimes]-t_start)/dt); 
spikes_ind_t(find(spikes_ind_t==0)) = 1;
%Spike Indices - cell
spikes_ind_c = [cells.index4spikes];

%Recount in overlapping bins for "boxcar" style (needed for MUAhist)
if exist('overlap','var')
    spikes_ind_c = repmat(spikes_ind_c,1,overlap);
    spikes_ind_t_temp = [];
    for o = 1:overlap
        spiketimeoffset = o-ceil(overlap/2);
        spikes_ind_t_temp = [spikes_ind_t_temp spikes_ind_t+spiketimeoffset];
    end
    spikes_ind_t = [];
    spikes_ind_t = spikes_ind_t_temp;
    %Remove negative t or overhanging spikes
    spikes_ind_c(find(spikes_ind_t<=0 |spikes_ind_t>=numts)) = [];
    spikes_ind_t(find(spikes_ind_t<=0 |spikes_ind_t>=numts)) = [];
end

%Spike Indices - convert to linear index
spikes_ind = sub2ind(size(spkmat), spikes_ind_t,spikes_ind_c);

%Add Spikes to bins
while spikes_ind
    spkmat(spikes_ind) = spkmat(spikes_ind)+1;
    %Remove full bins
    [uniquespikes, unisp_ind] = unique(spikes_ind);
    spikes_ind(unisp_ind) = [];
end

t = [0:size(spkmat,1)-1]'*dt+0.5*dt+t_start; %time vector (midpoint)

%spindices = [[cells.spiketimes]',[cells.index4spikes]'];

spikemat.data = spkmat;
spikemat.timestamps = t;
spikemat.dt = dt;

end

