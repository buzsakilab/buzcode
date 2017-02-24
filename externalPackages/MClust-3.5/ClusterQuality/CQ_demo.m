% load demo data
load('Data-R037-2003-09-04-TT07-CorticalCells.mat','T','WV','ClusterSpikes')

% create feature space 
Fet = Create_FeatureSpace(WV);

% Show clusters and quality
figure('Name','Cluster Quality Demo')
plot(Fet(:,1),Fet(:,2),'.k','MarkerSize',1);
xlabel('Energy on channel 1'); ylabel('Energy on channel 2');
hold on;

h = plot(NaN,NaN,'.b');

nC = length(ClusterSpikes);  % number of clusters
C_colors = {'r','b','g','m','c'};

disp('Press ENTER to step through the clusters...')
for iC = 1:nC
    % show the current cluster
    set(h,'Xdata',Fet(ClusterSpikes{iC},1),'Ydata',Fet(ClusterSpikes{iC},2),'Color',C_colors{iC});
    
    % calculate cluster quality
    CluSep = Cluster_Quality(Fet,ClusterSpikes{iC});
    % display the quality values
    CQ_string = sprintf('Cluster %1i: L_{ratio} = %5.4f  Isolation Distance = %4.1f',iC,CluSep.Lratio,CluSep.IsolationDist);
    title(CQ_string);
    
    CQ_string = sprintf('L_{ratio} = %5.4f  ID = %4.1f',CluSep.Lratio,CluSep.IsolationDist);
    text(1000,800 - 100*(iC),CQ_string,'Color',C_colors{iC},'FontSize',12)
    pause
    
    % leave the cluster on the plot shown in smaller points
    plot(Fet(ClusterSpikes{iC},1),Fet(ClusterSpikes{iC},2),'*','MarkerSize',1,'Color',C_colors{iC});
end
set(h,'Xdata',NaN,'Ydata',NaN);
title(' ');