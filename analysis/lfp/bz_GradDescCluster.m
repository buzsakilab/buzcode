function [ cluass,cluNRG ] = bz_GradDescCluster( cohmat, varargin )
%[cluass] = GradDescCluster(cohmat) clusters recording sites given a pairwise 
%similarity matrix using gradient descent to find minimize within-cluster 
%interaction energy. (i.e. maximize within-cluster coherence, see Berenyi 
%et al 2014 for details). Should be repeated a few times to confirm.
% 
%INPUT
%   cohmat      undirected pairwise connectivity matrix (coherence or other)
%   (optional parameters)
%   'numsteps'  number of steps (default: 5e5)
%   'numinit'   number of initial clusters, should be larger than the 
%               number of clusters you expect (default: 20)
%   'showplot'  true or false, to show the update plot as it goes
%   'brainmap'  cell array with two cells: {1} = linear indices, {2} = size map; 2 dim vector x and y dimensions of map
%   
%
%OUTPUT
%   cluass      cluster assignments for each 
%   cluNRG      energy (-mean within-cluster coherence) for each cluster
%
%From Berenyi et al. (2014). Large-scale, high-density (up to 512 
%   channels) recording of local circuits in behaving animals. Journal of 
%   Neurophysiology, 111(5), 1132?1149. http://doi.org/10.1152/jn.00785.2013
%
%Code implementation by DLevenstein 2016
%TO DO
%-Automated removal of small/scattered clusters if desired?
%% Parse the input parameters
parms = inputParser;
addParameter(parms,'numsteps',500000,@isnumeric);
addParameter(parms,'numinit',20,@isnumeric);
addParameter(parms,'showplot',true,@islogical);
addParameter(parms,'brainmap',{},@iscell);

parse(parms,varargin{:})
numsteps = parms.Results.numsteps;
numinit = parms.Results.numinit;
SHOWPLOT = parms.Results.showplot;
LinearInds = parms.Results.brainmap{1};
SizeVid = parms.Results.brainmap{2};

%% Initiate and Run
numsites = size(cohmat,1);
rng('shuffle') %Shuffle the random number generator seed... just in case
%Start from random initial assignments
cluass = randi(numinit,numsites,1);

if SHOWPLOT; figure; switchtracker = zeros(numsteps,1); end
for ss = 1:numsteps
    %The time ticker...
    if mod(ss,1000)==1
        display(['Step: ',num2str(ss),' of ',num2str(numsteps)])
        [~,clusort] = sort(cluass);
        map = ExpandVid(cluass, LinearInds, SizeVid);
        
        if SHOWPLOT
            subplot(2,2,1)
                imagesc(cohmat(clusort,clusort))
            subplot(2,2,2)
                plot(smooth(switchtracker(1:ss),500,'moving'))
                xlabel('Step #');ylabel('P(stay)')
                ylim([0 1])
            if ~isempty(LinearInds)
                subplot(2,2,3)
                imagesc(map)
            end;
            subplot(2,2,4)
                plot(ss,length(unique(cluass)),'o')
                hold on
                xlabel('Step #');ylabel('# Clusters')
            drawnow
            %hold on
        end
    end
    
    %Descend one step!
    cluass = DescOneStep(cluass,cohmat);
end

%Here: calculate Energy for each final cluster
finalclus = unique(cluass);
numfinalclus = length(finalclus);
for ff = 1:numfinalclus
    N = sum(cluass == finalclus(ff)); 
    %note: site i included in A and B
    E(ff) = -(1./N).*(sum(cohmat(cluass == finalclus(ff),cluass == finalclus(ff))));
end
cluNRG = [finalclus,E];


%% FUNCTION: OneStep
    function [newcluass] = DescOneStep(oldcluass,cohmat)
        clus = unique(oldcluass);
        numclus = length(clus);

        %Pick a random site and calculate the energy for it's current cluster
        sitepick = randi(length(oldcluass),1);
        clu_A = oldcluass(sitepick); %Cluster the current site is in       
        othersites = setdiff(1:length(oldcluass),sitepick);
        N_A = sum(oldcluass == clu_A)-1; %ALl OTHER elements of current cluster
        E_A = -(1./N_A).*sum(cohmat(sitepick,oldcluass(othersites)==clu_A));
        
        %If only one site in cluster, energy is 0
        if N_A == 0
            E_A = 0;
        end
        
        %Calculate energy gap for switching site to other clusters
        energygap = zeros(size(clus));
        for cc = 1:numclus
            clu_B = clus(cc);
            if clu_B==clu_A
                continue
            end
            N_B = sum(oldcluass == clu_B); 
            E_B = -(1./N_B).*(sum(cohmat(sitepick,oldcluass==clu_B)));
            energygap(cc) = E_A - E_B;
        end

        %Switch site to cluster with largest energy gap
        newcluass = oldcluass;
        newcluass(sitepick) = clus(energygap==max(energygap));
        
        %Keep Track of switching - 1 if the channel stays in its cluster
        if SHOWPLOT
            switchtracker(ss) = clus(energygap==max(energygap))==clu_A;
        end
    end
end

%% FUNCTION: ExpandVid

function map = ExpandVid(cluass, LinearInds, SizeVid)

%Intialize
map = single(nan(SizeVid(1)*SizeVid(2),1));

% Put tmpR values int initialized 2D matrix
map(LinearInds) = cluass;

% Reshape into video
map = reshape(map,SizeVid(1),SizeVid(2));
end
