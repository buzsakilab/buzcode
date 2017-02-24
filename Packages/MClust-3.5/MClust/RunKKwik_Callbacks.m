function FBP_Callbacks

% CBP_Callbacks
%
% Add a limit on the best projection of two clusters.
% 
% If cluster 0 is selected as the second cluster, then the program attempts
% to separate the cluster from the noise.  Takes a long time if you have a
% lot of spikes.
%
% INPUTS
% 
% NONE
%
% OUTPUTS
%
% NONE
% TO USE WITH MCLUST, put this in the MClust folder
% Also requires Find_Best_Projection.m

% NCST 2003
%

maxPossibleClusters = 20;

global MClust_CutOnBestProj
global MClust_FeatureTimestamps MClust_FeatureNames MClust_FeatureSources MClust_FDfn MClust_Clusters 
global MClust_TTdn MClust_TTfn MClust_TText MClust_ChannelValidity

newFeatures = get(findobj('Parent',findobj('Name','Run KlustaKwik'),'Tag','FeaturesUseListbox','TooltipString', ...
    'These features will be used for cluster separation.'),'String');

iClust = get(findobj('Tag','RunKlustaKwik'), 'UserData');
minClusters = str2double(get(findobj('Tag','RunKlustaKwikMinClust'),'String'));
maxClusters = str2double(get(findobj('Tag','RunKlustaKwikMaxClust'),'String'));

maxPossibleClusters = max([maxClusters maxPossibleClusters]);

figure(findobj('Tag','RunKlustaKwik'))
close;

if isempty(newFeatures)
    errordlg('No features chosen', 'RunKKwik Error', 'modal');
    return;
end

[f1 MClust_Clusters{iClust}] = FindInCluster(MClust_Clusters{iClust});

nFeatures = length(newFeatures);
nChannels = length(find(MClust_ChannelValidity));
FeatureData = repmat(nan,length(f1),nFeatures*nChannels);

[fpath fname fext] = fileparts(MClust_FDfn);

if isempty(fpath)
    fpath = pwd;
end 

for iX = 1:nFeatures
    FeatureToGet = newFeatures{iX};
    FDfn = fullfile(fpath, [fname '_' FeatureToGet '.fd']);
    if ~exist(FDfn)
        Get_MClust_FD_Defaults
        
        Write_fd_file(fpath, fullfile(MClust_TTdn,[MClust_TTfn MClust_TText]), ...
            newFeatures, MClust_ChannelValidity, record_block_size, template_matching, NormalizeFDYN)
    end
    temp = load(FDfn,'-mat');
	f_CH = find(ismember(find(MClust_ChannelValidity),find(temp.ChannelValidity)));
    FeatureData(:,(1:nChannels) + nChannels*(iX - 1)) = temp.FeatureData(f1,f_CH);
    DisplayProgress(iX,nFeatures);
end

FETfn = fullfile(fpath,[fname '-Cluster' num2str(iClust)]);
WriteFeatureData2TextFile(FETfn,FeatureData);

KlustaKwikPath = which('KlustaKwik.exe');
if ~isempty(KlustaKwikPath)
    disp(['Using ' KlustaKwikPath]);
else
    disp('Did not find KlustaKwik.exe.');
    return
end     

file_no = 1;
parameter_string = ['-MinClusters ' num2str(minClusters) ...
        ' -MaxClusters ' num2str(maxClusters) ...
        ' -MaxPossibleClusters ' num2str(maxPossibleClusters) ' -UseFeatures ' num2str(repmat(1,size(FeatureData,2),1))'];
COMMAND = ['! ' KlustaKwikPath ' "' FETfn '" ' num2str(file_no) ' ' parameter_string ];
disp(FETfn);
disp(['Number of spikes: ' num2str(size(FeatureData,1))]);
%                disp(['Estimated time required to run KlustaKwik.exe on this file is ' num2str(EstDuration(i)*60) ' seconds (or ' num2str(EstDuration(i)/60) 'hours)']);
COMD_output = evalc(COMMAND);
% Find the output that you are going to display
COMD_filebase = findstr(COMD_output,'FileBase');
COMD_dim = findstr(COMD_output,'dimension');
COMD_time = findstr(COMD_output,'That took');

try
    disp(COMD_output([COMD_filebase:COMD_dim + 13,COMD_time:end]));
catch
    disp(['Did not find output that was searched for, all output shown:  ' COMD_output]);
end

% diary on
% disp(COMD_output)
% diary off

disp(' ')

clusters = load([FETfn '.clu.1']);
clusters = clusters(2:end);
unique_clusters = unique(clusters);
n_clusters = length(unique_clusters);

for iC = 1:n_clusters
    idx = find(clusters == unique_clusters(iC));
	MClust_Clusters{end + 1} = precut(['Cluster ' num2str(iClust) '_sub' num2str(iC)]);
	MClust_Clusters{end} = AddIndices(MClust_Clusters{end}, f1(idx));
end

figHandle = findobj('Type','figure','Tag', 'ClusterCutWindow');
MClustCutterClearClusterKeys(figHandle);
MClustCutterRedrawClusterKeys(figHandle, max(0,length(MClust_Clusters)-16));
MClustCutterCallbacks('RedrawAxes')

% ClusterCutWindow = findobj('Type','figure','Tag', 'ClusterCutWindow');
% if ~isempty(ClusterCutWindow)
%     close(ClusterCutWindow);
% end
% 
% CHDrawingAxisWindow = findobj('Type','figure','Tag', 'CHDrawingAxisWindow');
% if ~isempty(CHDrawingAxisWindow)
%     close(CHDrawingAxisWindow);
% end

fn = FindFiles([FETfn '*']);
for iFN = 1:length(fn)
    eval(['! del ' fn{iFN}]);
end


%===============================================================================
function WriteFeatureData2TextFile(file_name, FeatureData)
%
% write featuredata from memory to a text file for input into KlustaKwick.exe
%
file_no = 1;
fid = fopen([ file_name '.fet.' num2str(file_no)],'w');
[n_points, n_features] = size(FeatureData);
fprintf(fid,'%3d \n',n_features);
for ii = 1:n_points
    fprintf(fid,'%f\t',FeatureData(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
