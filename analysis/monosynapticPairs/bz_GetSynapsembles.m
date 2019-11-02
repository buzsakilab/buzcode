function bz_GetSynapsembles(varargin)
% bz_GetSynapsembles - identify E->I synaptic assemblies
%
% USAGE
%
%    bz_GetSynapsembles(varargin)
%
% INPUTS
%
%    basepath               -path to recording (where .dat/.clu/etc files are)
%    cellTypeSource         -file format for loading cell types
%                           (default='SM' ; allowed formats {'buzcode', 'PP', 'SM'})
%    fitCoupling            -logical (default=false) to fit the long term
%                           coupling; if coupling has already been fit, this will
%                           refit and overwrite existing files
%    saveMat                -logical (default=false) to save in buzcode format
%    forcePlasticity        -logical (default=false) run bz_GetSynapsembles.m for 
%                            extracting fluctuations in spike transmission only 
%                            (up to bz_fitPoissPlasticity.m)
%
% OUTPUTS
%
% 1) Find/load connected pairs - find significant high frequency, delayed
%      synchrony (English, McKenzie et al., Neuron  2017)
% 2) Compute/load the temporal fluctuations in synaptic coupling for each pair
%      while  regressing out slow changes in post-synaptic rates
% 3) Identify synapses that covary together - synapsembles
%      PCA then ICA (method described in Lopes-dos-Santos et al., 2013)
%
% NOTES 
%   - the function is written to support file formats in datasets that have
%     looked at monosynaptic connections
%   - for now, the function allows non-buzcode formats for some of the
%     necessary inputs (e.g., cell types)
% 
% Roman Huszar & Sam McKenzie, August 2019

%% STEP 0 PROCESS INPUTS

% Parse / save user inputs
p = inputParser;

addParameter(p, 'cellTypeSource', 'buzcode', @ischar)
addParameter(p, 'basepath',pwd,@isstr);
addParameter(p, 'fitCoupling',false,@islogical);
addParameter(p, 'plotFlag',false,@islogical); 
addParameter(p, 'forcePlasticity',false,@islogical)

parse(p,varargin{:})

cellTypeSource = p.Results.cellTypeSource;
basepath = p.Results.basepath;
fitCoupling = p.Results.fitCoupling;
plotFlag        = p.Results.plotFlag;
forcePlasticity = p.Results.forcePlasticity;

% If 'monosyn_plast' folder does not exist or is empty, fit long term coupling
glmfits_datadir = fullfile(basepath,'monosyn_plast');
if ~exist(glmfits_datadir, 'dir') || ( exist(glmfits_datadir, 'dir') && 2 == length(dir(glmfits_datadir)) )
    fitCoupling = true;
end

basename = bz_BasenameFromBasepath(basepath);

% Get cell types
if strcmp(cellTypeSource, 'buzcode')   % buzcode celltype formatting
    try
        load(fullfile(basepath,[basename '.CellClass.cellinfo.mat']))
    catch
        error('basename.CellClass.cellinfo.mat could not be found in basepath...')
    end
    pyr = find( CellClass.pE );
    int = find( CellClass.pI );
elseif strcmp(cellTypeSource, 'SM')    % Sam McKenzie's celltype formatting
    try
        v = load(fullfile(basepath,[basename '_CellParams.mat']));
    catch
        error('basename_CellParams.mat could not be found in basepath...')
    end
    pyr = find( [v.CellParams.cellType] == 1 );
    int = find( [v.CellParams.cellType] > 1 );
elseif strcmp(cellTypeSource, 'PP')    % Peter Petersen's celltype formatting
    try
        load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    catch
        error('basename.cell_metrics.cellinfo.mat could not be found in basepath...')
    end
    pyr = find( cellfun( @(x) strcmp(x, 'Pyramidal Cell'), cell_metrics.putativeCellType ) );
    int = find( cellfun( @(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType ) );
else
    error('Incorrect value for property ''cellTypeSource'' ');
end

%% STEP 1 FIND CONNECTED PAIRS (PYR-->INT)
%  English, McKenzie et al., 2017 

if fitCoupling

% Find monosynaptically connected pairs
mono_path = [basepath filesep basename '.mono_res.cellinfo.mat'];
if strcmp(cellTypeSource, 'SM')
    mono_res = v.mono_res;
elseif exist(mono_path, 'file')
    load(mono_path);
else
    fprintf('Detecting monosynaptic connections...\n')
    mono_res = bz_GetMonoSynapticallyConnected(basepath,'plot',false, 'saveMat', false);
end

% Only take connections between pyr --> int, other significant hits are
% either high freq. interneuron synchrony likely due to common input or
% spike sorting errors where spikes that are part of a burst a shared
% across units
if ~isfield(mono_res, 'pyr2int')
    kp = ismember(mono_res.sig_con(:,1),pyr) & ismember(mono_res.sig_con(:,2),int);
    mono_res.pyr2int = mono_res.sig_con(kp,:);
    % Store mono_res.cellinfo file that includes 
    save( fullfile(basepath, [basename '.mono_res.cellinfo.mat']), 'mono_res')
end

end

%% 2) ESTIMATE FLUCTUATIONS IN EXCESS SYNCHRONY OVER TIME

if fitCoupling

% Fit Poisson GLM to account for long term coupling
for ii = 1:size(mono_res.pyr2int,1)
    
    bz_fitPoissPlasticity(basepath, mono_res.pyr2int(ii,1), mono_res.pyr2int(ii,2))
    
end

end

% Return if we only want to extract time-varying spike transmission
if forcePlasticity
   return 
end

%% 3) DETECT SYNAPSEMBLES

% Load all data from step 2
fils= dir(glmfits_datadir);
fils = {fils.name};
kp = cellfun(@any,regexp(fils',[basename '_p']));
fils = fils(kp);

% Load all spike transmission timeseries
load(fullfile(glmfits_datadir,fils{1}), 'ts_lastelement', 'down_dt','dt')
ts = 0:dt:ts_lastelement; ts_down = ts(1:round(down_dt/dt):end);
spktrans_t = nan(length(fils),  length(ts_down));
pairID = nan(length(fils), 2);

fprintf('Loading spike transmission timeseries... \n')
for ii = 1:length(fils)
    
    v = load(fullfile(glmfits_datadir,fils{ii}), 'ltc_t', 'preID', 'postID');
    spktrans_t(ii,:) = v.ltc_t;
    pairID(ii,:) = [v.preID v.postID];
    
end

% Get rid of pairs with bad fits (nans / infs)
pairID( any(isnan(spktrans_t),2) | any(isinf(spktrans_t),2) , : ) = [];
spktrans_t( any(isnan(spktrans_t),2) | any(isinf(spktrans_t),2) , : ) = [];

% Parameters for assembly detection
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;

fprintf('Detecting synapsembles (PCA/ICA)... \n')
assemblyTemplates = assembly_patterns(spktrans_t,opts);

% Project each datapoint in pre / post stim epoch onto each assembly
% template to get assembly expression strength over time
assemblyExpressions = assembly_activity(assemblyTemplates, spktrans_t);

% Save output of synapsemble analysis
fprintf('Saving synapsembles... \n')
save( fullfile( glmfits_datadir, 'synapsembles.mat' ), 'assemblyTemplates', 'assemblyExpressions', '-v7.3')

%% 4) PLOT

if plotFlag

close all

figure
set(gcf,'Position',[ 1000 833 1109 505 ])

% Find sorting that maximally decorelates synapsemble weigths
D = pdist(assemblyTemplates,'correlation');
Z = linkage(D,'average');
leafOrder = optimalleaforder(Z,D);

% Find which pairs best correlate with which synapsemble
wc = corr(zscore([assemblyExpressions;spktrans_t(leafOrder,:)])');
nAssemblies = size(assemblyExpressions,1);
wc = wc(nAssemblies+1:end,1:nAssemblies);
[~,b] = max(wc);
[~,b]  = sort(b);
aw = assemblyExpressions(b,:);

% Plot spike transmission timeseries
subplot(2,1,1)
imagesc(ts_down / 60,[],zscore(spktrans_t(leafOrder,:),[],2),[-5 5])
colorbar
xlabel('Time (min)', 'fontsize', 14)
ylabel('PYR-INT pair', 'fontsize', 14)

% Plot synapsemble expression strength timeseries
subplot(2,1,2)
imagesc(ts_down / 60,[],zscore(aw,[],2),[-5 5])
xlabel('Time (min)', 'fontsize', 14)
ylabel('Synapsemble', 'fontsize', 14)
colormap(gca,bluewhitered)
colorbar

set(gca,'ydir','normal')

end

%%

end
