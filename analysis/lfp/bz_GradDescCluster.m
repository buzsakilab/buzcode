function [ cluass,cluNRG ] = bz_GradDescCluster( cohmat, varargin )
%[cluass] = GradDescCluster(cohmat) clusters recording sites given a pairwise 
%similarity matrix using gradient descent to find minimize within-cluster 
%interaction energy. (i.e. maximize within-cluster coherence, see Berenyi 
%et al 2014 for details). Should be repeated a few times to confirm.
% 
%INPUT
%   cohmat      undirected pairwise connectivity matrix (coherence or other)
%   (optional parameters)
%   'numsteps'      max. number of steps (default: 5e5)
%   'numinit'       number of initial clusters, should be larger than the 
%                   number of clusters you expect (default: 20)
%   'showplot'      true or false, to show the update plot as it goes
%   'stopthresh'    stability threshold to stop descending. probability of
%                   selected site changing cluster identity, P(switch)
%                   (default: 0.001)
%   'stopwin'       window of time for which the P(switch) must stay below
%                   stopthresh (default: 10000 steps)
%
%OUTPUT
%   cluass      cluster assignments for each 
%   cluNRG      neg. energy (i.e. mean coherence) for each cluster
%
%From Berenyi et al. (2014). Large-scale, high-density (up to 512 
%   channels) recording of local circuits in behaving animals. Journal of 
%   Neurophysiology, 111(5), 1132?1149. http://doi.org/10.1152/jn.00785.2013
%
%Code implementation by DLevenstein 2016
%% Parse the input parameters
parms = inputParser;
addParameter(parms,'numsteps',500000,@isnumeric);
addParameter(parms,'numinit',20,@isnumeric);
addParameter(parms,'showplot',true,@islogical);
addParameter(parms,'stopthresh',0.001,@(x) x>0 && x<=1)
addParameter(parms,'stopwin',10000,@isnumeric)

parse(parms,varargin{:})
numsteps = parms.Results.numsteps;
numinit = parms.Results.numinit;
SHOWPLOT = parms.Results.showplot;
stopthresh = parms.Results.stopthresh;
stopwin = parms.Results.stopwin;

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
        
        if SHOWPLOT
            subplot(2,2,1)
                imagesc(cohmat(clusort,clusort))
            subplot(2,2,2)
                plot(log10(smooth(switchtracker(1:ss),1000,'moving')))
                hold on
                plot(get(gca,'xlim'),log10(stopthresh).*[1 1],'r--')
                xlabel('Step #');ylabel('P(switch)')
                ylim([log10(stopthresh)-1 0])
                LogScale('y',10)
                hold off
            subplot(2,2,4)
                plot(ss,length(unique(cluass)),'ko')
                hold on
                xlabel('Step #');ylabel('# Clusters')
            drawnow
            %hold on
        end
    end
    
    %Descend one step!
    [cluass,switchtracker(ss)] = DescOneStep(cluass,cohmat);
    
    %Decide to exit the loop - if you've been consistent for a certain
    %window
    if mean(switchtracker(max(ss-stopwin,1):ss))<=stopthresh
        display('Clustering has descended to the bottom of the (local) pit!')
        break
    end
end

%Here: calculate -Energy (mean coherence)for each final clusters
finalclus = unique(cluass);
numfinalclus = length(finalclus);
for ff = 1:numfinalclus
    %N = sum(cluass == finalclus(ff)); 
    clusmat = cohmat(cluass == finalclus(ff),cluass == finalclus(ff));
    cluspairs = triu(clusmat,1);
    E(ff) = mean(cluspairs(:));
end
cluNRG = [finalclus,E'];


%% FUNCTION: OneStep
    function [newcluass,diditswitch] = DescOneStep(oldcluass,cohmat)
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
        diditswitch = clus(energygap==max(energygap))~=clu_A;
    end
end

