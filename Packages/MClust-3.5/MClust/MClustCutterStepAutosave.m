function MClustCutterStepAutosave(force)

%========================================
% MClustCutterStepAutosave(force)
% steps autosave using global variable MClust_autosave
% if force is true then autosaves anyway

if nargin == 0; force = 0; end

global MClust_Colors MClust_TTfn MClust_Clusters MClust_FeatureNames MClust_ChannelValidity
global  featureindex featuresToUse
global MClust_Autosave

MClust_Autosave = MClust_Autosave-1;
autosaveCounterHandle = findobj('Tag', 'Autosave');

if MClust_Autosave == 0 || force
    set(autosaveCounterHandle, 'String', 'Autosaving ...');
    % get components
    [basedn, basefn, ext] = fileparts(MClust_TTfn);
    featureToUseHandle = findobj('Tag', 'FeaturesUseListbox');
    featuresToUse = get(featureToUseHandle, 'String');
    % save defaults
    save(fullfile(basedn,'autodflts.mclust'), 'MClust_Colors', 'featuresToUse', '-mat');
    % save clusters
    save(fullfile(basedn,'autosave.clusters'), 'MClust_Clusters', 'MClust_FeatureNames', ...
         'MClust_ChannelValidity','MClust_Colors','featuresToUse','featureindex','-mat');
    % reset counter
    MClust_Autosave = 10;
    set(autosaveCounterHandle, 'String', 'Autosaved');
end

set(autosaveCounterHandle, 'String', ['Autosave in ' num2str(MClust_Autosave)]);