function [ cluass,cluNRG ] = bz_GradDescCluster( simMat, varargin )
%[cluass] = bz_GradDescCluster(simMat) clusters recording sites given a pairwise 
%similarity matrix using gradient descent to find minimize within-cluster 
%interaction energy. (i.e. maximize within-cluster coherence, see Berenyi 
%et al 2014 for details). 
% 
%INPUT
%   simMat          undirected similarity matrix (coherence, correlation, or other)
%
%   (optional parameters)
%   'numsteps'      number of steps (default: 5e5)
%   'numinit'       number of initial clusters, should be larger than the 
%                   number of clusters you expect (default: 20)
%   'stopthresh'    stability threshold to stop descending. probability of
%                   selected site changing cluster identity, P(switch)
%                   (default: 0.001)
%   'stopwin'       window of time for which the P(switch) must stay below
%                   stopthresh (default: 10000 steps)   
%   'showplot'      true or false, to show the update plot as it goes
%   'brainmap'      cell array with two cells: {1} = linear indices, {2} = size map; 2 dim vector x and y dimensions of map

%OUTPUT
%   cluass      cluster assignments for each 
%   cluNRG      negative energy (i.e. mean within-cluster coherence) for each cluster
%
%From Berenyi et al. (2014). Large-scale, high-density (up to 512 
%   channels) recording of local circuits in behaving animals. Journal of 
%   Neurophysiology, 111(5), 1132?1149. http://doi.org/10.1152/jn.00785.2013
%
%Code implementation by DLevenstein 2016. Updates DL and RS May 2017.
%
% TO DO
% - make brainmap general to any 2D ephys/imaging data 
% - develop method for stopping num iterations
% - add consensus clustering option? 

%% Parse the input parameters
parms = inputParser;
addParameter(parms,'numsteps',500000,@isnumeric);
addParameter(parms,'numinit',20,@isnumeric);
addParameter(parms,'showplot',true,@islogical);
addParameter(parms,'brainmap',{},@iscell); 
addParameter(parms,'stopthresh',0.001,@(x) x>0 && x<=1)
addParameter(parms,'stopwin',10000,@isnumeric)

parse(parms,varargin{:})
numsteps = parms.Results.numsteps;
numinit = parms.Results.numinit;
SHOWPLOT = parms.Results.showplot;
LinearInds = parms.Results.brainmap{1};
SizeVid = parms.Results.brainmap{2};
stopthresh = parms.Results.stopthresh;
stopwin = parms.Results.stopwin;

%% Initialize and Run
numsites = size(simMat,1);
rng('shuffle') %Shuffle the random number generator seed... just in case
cluass = randi(numinit,numsites,1); %Start from random initial assignments

if SHOWPLOT; figure; switchtracker = zeros(numsteps,1); end 
for ss = 1:numsteps
    %The time ticker...
    if mod(ss,1000)==1
        display(['Step: ',num2str(ss),' of ',num2str(numsteps)])
        [~,clusort] = sort(cluass);
        
        % Visuals
        if SHOWPLOT
            % plot changing clusters over time
            subplot(2,2,1) 
                imagesc(simMat(clusort,clusort)) 
            %plot switching behavior
            subplot(2,2,2)
                plot(log10(smooth(switchtracker(1:ss),1000,'moving'))) 
                hold on
                plot(get(gca,'xlim'),log10(stopthresh).*[1 1],'r--')
                xlabel('Step #');ylabel('P(switch)')
                ylim([log10(stopthresh)-1 0])
                LogScale('y',10)
                hold off
            % plot map    
            if ~isempty(LinearInds) 
                map = ExpandVid(cluass, LinearInds, SizeVid);
                subplot(2,2,3)
                imagesc(map)
            end;
            %plot number of clusters
            subplot(2,2,4)
                plot(ss,length(unique(cluass)),'ko') %number clusters
                hold on
                xlabel('Step #');ylabel('# Clusters')
                drawnow
        end;
    end;
    
    % Fx meat
    [cluass,switchtracker(ss)] = DescOneStep(cluass,simMat);
    
    %Decide to exit the loop - if you've been consistent for a certain window
    if mean(switchtracker(max(ss-stopwin,1):ss)) <= stopthresh
        disp('Clustering has descended to the bottom of the (local) pit!')
        break
    end
end

%Calculate -Energy (mean similarity) for each final cluster
finalclus = unique(cluass);
numfinalclus = length(finalclus);
for ff = 1:numfinalclus
    %N = sum(cluass == finalclus(ff)); 
    clusmat = simMat(cluass == finalclus(ff),cluass == finalclus(ff));
    cluspairs = triu(clusmat,1);
    E(ff) = mean(cluspairs(:));
end
cluNRG = [finalclus,E'];



%% FUNCTION: OneStep
    function [newcluass,diditswitch] = DescOneStep(oldcluass,simMat)
        clus = unique(oldcluass);
        numclus = length(clus);

        %Pick a random site and calculate the energy for it's current cluster
        sitepick = randi(length(oldcluass),1);
        clu_A = oldcluass(sitepick); %Cluster the current site is in       
        othersites = setdiff(1:length(oldcluass),sitepick);
        N_A = sum(oldcluass == clu_A)-1; %ALl OTHER elements of current cluster
        E_A = -(1./N_A).*sum(simMat(sitepick,oldcluass(othersites)==clu_A));
        
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
            E_B = -(1./N_B).*(sum(simMat(sitepick,oldcluass==clu_B)));
            energygap(cc) = E_A - E_B;
        end

        %Switch site to cluster with largest energy gap
        newcluass = oldcluass;
        newcluass(sitepick) = clus(energygap==max(energygap));
        
        %Keep Track of switching - 1 if the channel stays in its cluster
        diditswitch = clus(energygap==max(energygap))~=clu_A;
    end
end

%% FUNCTION: ExpandVid

function map = ExpandVid(cluass, LinearInds, SizeVid)

    % Intialize
    map = single(nan(SizeVid(1)*SizeVid(2),1));

    % Put tmpR values int initialized 2D matrix
    map(LinearInds) = cluass;

    % Reshape into video
    map = reshape(map,SizeVid(1),SizeVid(2));

end
