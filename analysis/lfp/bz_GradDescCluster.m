function [ cluass ] = GradDescCluster( cohmat, varargin )
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
%
%OUTPUT
%   cluass      cluster assignments for each 
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

parse(parms,varargin{:})
numsteps = parms.Results.numsteps;
numinit = parms.Results.numinit;

%% Initiate and Run
numsites = size(cohmat,1);
%Start from random initial assignments
cluass = randi(numinit,numsites,1);

for ss = 1:numsteps
    cluass = DescOneStep(cluass,cohmat);
end


%% FUNCTION: OneStep
    function [newcluass] = DescOneStep(oldcluass,cohmat)
        clus = unique(oldcluass);
        numclus = length(clus);

        %Pick a random site and calculate the energy for it's current cluster
        sitepick = randi(length(oldcluass),1);
        %othersites = setdiff(1:length(oldcluass),sitepick);
        clu_A = oldcluass(sitepick);
        N_A = sum(oldcluass == clu_A);
        E_A = -(1./N_A).*sum(cohmat(sitepick,oldcluass==clu_A));

        %Calculate energy gap for switching site to other clusters
        energygap = zeros(size(clus));
        for cc = 1:numclus
            clu_B = clus(cc);
            if clu_B==clu_A
                continue
            end
            N_B = sum(oldcluass == clu_B)+1; 
            %note: site i included in A and B
            E_B = -(1./N_B).*(sum(cohmat(sitepick,[find(oldcluass==clu_B);sitepick])));
            energygap(cc) = E_A - E_B;
        end

        %Switch site to cluster with largest energy gap
        newcluass = oldcluass;
        newcluass(sitepick) = clus(energygap==max(energygap));

    end
end

