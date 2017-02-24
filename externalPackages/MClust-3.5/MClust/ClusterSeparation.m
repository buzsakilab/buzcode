function [L_Extra,L_Ratio,IsolationDist,Dists,df,L_Intra,C] = ClusterSeparation(f,filename,ChannelValidity,iClust,varargin)

% [L_Extra,L_Ratio,IsolationDist,Dists,df] = ClusterSeparation(f,filename,ChannelValidity,iClust);
%
% ClusterSeparation.m
%
% Evaluates the separation of a cluster from other spikes/noise on the same
% tetrode/electrode.  Returns 4 values related to cluster quality, and the
% pdfs of distances in the cluster and those outside of the cluster.  If no
% outputs are given (nargout == 0), creates a figure showing the
% separation.
% 
% INPUTS: f -- index of points in the cluster
%         filename -- name of data file from which f was created
%         ChannelValidity -- 4 element vector specifying which channels to
%           use.
%         iClust -- Cluster number (just for titling output figure)
% TO USE WITH MCLUST, put this in the MClust/ClusterOptions folder
% ncst 14 Nov 02 with strong input from adr and jjackson
% STATUS: beta (in use)
% adr 17 aug 03 skip with warning if toolbox unavail

% warning off MATLAB:divideByZero

global MClust_ClusterSeparationFeatures

if ~license('test', 'Statistics_Toolbox') || isempty(which('mahal'))
	 % needs mahal from Stats toolbox
	warning('MCLUST:ToolboxUnavailable', 'Skipping cluster separation. Stats toolbox unavailable.');
	L_Extra = nan;
	L_Ratio = nan;
	IsolationDist = nan;
	Dists = [];
	df = nan;
	L_Intra = nan;
	C = nan;
	return
end

MaxRecords = 400000;     % maximum number of records to load at one time
Dists = [];              % Mahalanobis distance of each spike to the center of points defined by f
df = [];                 % degrees of freedom for chi^2
TTFileName = [];         % name of original tetrode file 
doplot = 0;              % plot output
n_records = [];          % number of records in TTFileName
FeaturesToGet = {};      % Features to use in calculating separation
FD = [];                 % Feature data matrix
extract_varargin;

if nargin < 3
	ChannelValidity = [1 1 1 1];
end

if nargin < 4
	iClust = ' ';
end

NewDir = 0;

if isempty(FD)  % If no feature data matrix was supplied
	if isempty(Dists) || isempty(df)  % and no Mahalanobis distance was supplied
		[fpath fname fext] = fileparts(filename);	
		
		if isempty(TTFileName)
			TTFileName = [pwd filesep fname fext];
			if ~isempty(fpath)
				pushdir(fpath);
				TTFileName = [fpath filesep fname fext];
				NewDir = 1;
			end
		end
		
		if exist('FD','dir')
			FDdn = [pwd filesep 'FD'];
		else
			FDdn = pwd;
		end
		
		
		if ~isempty(MClust_ClusterSeparationFeatures)
			FeaturesToGet = MClust_ClusterSeparationFeatures;
		else
			disp('No MClust_ClusterSeparationFeatures variable found: using default features')
			FeaturesToGet = {'Energy'};
		end
		
		nFeatures = length(FeaturesToGet);
		
		FD = [];
		FDfiles = FindFiles([fname '_*.fd']);
		if ~isempty(FDfiles)
			temp = load(FDfiles{1},'-mat');
			if ~all(ChannelValidity == temp.ChannelValidity) && nargin < 3
				disp('Using Channel Validity from FD files')
				ChannelValidity = temp.ChannelValidity;
			end
		end
		
		% Find the size of the TT file
		if isempty(n_records)
			n_records = MClust_CountSpikes(TTFileName);
		end
		
		% Length(f) needs to be longer than the number of columns in FD
		if length(f) < length(find(ChannelValidity))*length(FeaturesToGet)
			disp('Not enough spikes to calculate mahalanobis distance')
% 			L_Extra,L_Ratio,IsolationDist,Dists,df,L_Intra,C
			Dists = repmat(nan,n_records,1);
			L_Ratio = NaN;
			L_Extra = nan;
			L_Intra = nan;
			IsolationDist = nan;
			df = length(find(ChannelValidity))*length(FeaturesToGet);
			pdfs = repmat(nan,100,3);
			pdfs(:,1) = (1:100)';
			return
		end
		
		CalcFeatures = 0;
		for iF = 1:nFeatures
			FDFileName = fullfile(FDdn, [fname '_' FeaturesToGet{iF} '.fd']);
			if ~exist(FDFileName, 'file') 
				CalcFeatures = 1;
			end
		end
		if CalcFeatures
			Get_MClust_FD_Defaults
			Write_fd_file(FDdn, TTFileName, ...
				FeaturesToGet, ChannelValidity, record_block_size, template_matching, NormalizeFDYN)
		end          
		
		if n_records <= MaxRecords
			for iF = 1:nFeatures
				FDFileName = fullfile(FDdn, [fname '_' FeaturesToGet{iF} '.fd']);
				temp = load(FDFileName,'-mat');
				if size(temp.FeatureData,2) == size(ChannelValidity,2)
					FD = [FD temp.FeatureData(:,find(ChannelValidity))];
				else
					FD = [FD temp.FeatureData];
				end
			end
			Dists = mahal(FD,FD(f,:)); 
		else
			Dists = [];
			FD_f = [];
			% Find feature data points for spikes in the cluster
			for iF = 1:nFeatures
				FDFileName = fullfile(FDdn, [fname '_' FeaturesToGet{iF} '.fd']);
% 				if ~exist(FDFileName) 
% 					record_block_size = 80000;  % maximum number of spikes to load into memory at one time
% 					template_matching = 0;
% 					NormalizeYN = 'yes';
% 					disp(' ')
% 					disp(['FD file ' FDFileName '.fd does not exist; creating...'])
% 					disp(' ')
% 					Write_fd_file(FDFileName, TTFileName, ...
% 						FeaturesToGet(iF), ChannelValidity, record_block_size, template_matching, NormalizeYN)
% 				end          
				temp = load(FDFileName,'-mat');
				if size(temp.FeatureData,2) == size(ChannelValidity,2)
					FD_f = [FD_f temp.FeatureData(f,ChannelValidity)];
				else
					FD_f = [FD_f temp.FeatureData(f,:)];
				end
			end
		
			nBlocks = ceil(n_records/MaxRecords);
			for iB = 1:nBlocks
				StartRecord = MaxRecords*(iB - 1) + 1;
				if iB ~= nBlocks
					EndRecord = MaxRecords*(iB);
				else
					EndRecord = n_records;
				end
				FD = [];
				for iF = 1:nFeatures
					FDFileName = fullfile(FDdn, [fname '_' FeaturesToGet{iF} '.fd']);
					temp = load(FDFileName,'-mat');
					if size(temp.FeatureData,2) == size(ChannelValidity,2)
						FD = [FD temp.FeatureData(StartRecord:EndRecord,ChannelValidity)];
					else
						FD = [FD temp.FeatureData(StartRecord:EndRecord,:)];
					end
				end
				Dists = [Dists; mahal(FD,FD_f)];
			end
		end
	else
		n_records = length(Dists);
		FD = repmat(nan,1,df);
	end
else
	Dists = mahal(FD,FD(f,:));
	n_records = length(Dists);
end

f_out = repmat(1,n_records,1);
f_out(f) = 0;
f_out = find(f_out);

IntraDist = Dists(f);
ExtraDist = Dists(f_out); 

%------------------------------------
% L_Ratio
% take distances from mahal.m to be distributed as chi^2 with df
% degrees of freedom [df = (# features)*(# channels)]
% L_Extra is the sum of the probability that spikes from
% outside the cluster (f_out) are in the probability density of the
% cluster (f).
if isempty(df)
	df = size(FD,2);
end

% Original L_Intra/L_Extra formulation
% L_Intra = sum(chi2pdf(Dists(f),df));
% L_Extra = sum(chi2pdf(Dists(f_out),df));

% Jadin suggested 1 - chi2pdf
L_Intra = sum(1-chi2cdf(Dists(f),df));
L_Extra = sum(1-chi2cdf(Dists(f_out),df));

L_Ratio = L_Extra/length(f);
%-------------------------------------------

%-------------------------------------------
% Ken Harris' cluster quality measure
IsolationDist = IsolationDistance(FD,f,Dists);
%-------------------------------------------

bins_In = logspace(log10(min(Dists)),log10(max(Dists)),floor(2*sqrt(length(IntraDist))));
H_In = hist(Dists(f),bins_In);
bins_Out = logspace(log10(min(Dists)),log10(max(Dists)),floor(2*sqrt(length(ExtraDist))));
H_Out = hist(Dists(f_out),bins_Out);

if nargout == 0 || doplot
	ScSize = get(0, 'ScreenSize');
	figure('Name',['Separation: Cluster ' num2str(iClust)],'Position',[ ScSize(3)/4 ScSize(4)*0.1...
			ScSize(3)/2 ScSize(4)*0.8]);
	
	%------------------------------------------
	% Text summary 
	subplot(2,2,1)	
	msgstr = {};
	msgstr{end+1} = ['Cluster ' num2str(iClust)];
	msgstr{end+1} = ' ';
	msgstr{end+1} = sprintf('L-Extra                = %5.3f',L_Extra);
	msgstr{end+1} = sprintf('L-Ratio                = %5.4f',L_Ratio);
% 	msgstr{end+1} = sprintf('Chi fit (corr)          = %3.2f',C);
	msgstr{end+1} = sprintf('Isolation Distance = %4.1f',IsolationDist);
	msgstr{end+1} = ' ';
	msgstr{end+1} = 'Features: ';
	for iF = 1:length(FeaturesToGet)
		msgstr{end+1} = ['   ' FeaturesToGet{iF}];
	end

	axis off;
	h=text(0,0.5, msgstr);
	
	%------------------------------------------
	% Histograms of distances inside (red) and outside (blue) of cluster
	subplot(2,2,2)
	hold on
	[AX,H1,H2] = plotyy(bins_In,H_In,bins_Out,H_Out); 
	axes(AX(1)); ylabel('In Cluster');
	axes(AX(2)); ylabel('Out of Cluster');
	
	set(H1,'LineWidth',2,'Color','r');
	set(H2,'LineWidth',2,'Color','b','LineStyle','-.');
	legend([H1 H2],'In Cluster','Out of Cluster')
	title('Raw Count')
	set(AX(1),'Xtick',logspace(0,floor(log10(max(Dists))),floor(log10(max(Dists)))+1),...
		'Xlim',[0 max(Dists)],'Xscale','log','Ycolor','r')
	set(AX(2),'Xtick',logspace(0,floor(log10(max(Dists))),floor(log10(max(Dists)))+1),...
		'Xlim',[0 max(Dists)],'Xscale','log','Ycolor','b')
	xlabel('Mahalanobis distance')
	
	%------------------------------------------
	% PDFs of distances inside (red) and outside (blue) of cluster as well
	% as the chi squared pdf (black) for points within the cluster.  If the
	% red line and black line are the same, then ChiVal is a valid
	% estimate.
	subplot(2,2,3)
	hold on
	mxDist = chi2inv(0.9999,df);
	bins = linspace(0,mxDist,round(1.5*sqrt(length(IntraDist))));
	H_In = hist(IntraDist(IntraDist < mxDist),bins);
	H_In = H_In./sum(H_In);
	
	plot(bins,H_In./sum(H_In),'r','LineWidth',2);
% 	plot(bins_Out,H_Out./sum(H_Out),'-.b','LineWidth',2);

	% Adding 2 to df seems to adjust for the use of log scale bins.
	plot(bins,chi2pdf(bins,df)/sum(chi2pdf(bins,df)), '--k','LineWidth',2)
	
% 	legend('In Cluster','Out of Cluster','Chi^2 fit')
	legend('In Cluster','Chi^2 fit')
	title('PDF')
	set(gca,'Xlim',[0 mxDist])
	xlabel('Mahalanobis distance')
	
	%------------------------------------------
	% BoxPlot of Distances inside and outside of the cluster
	subplot(2,2,4)
	boxplot([IntraDist; ExtraDist],[repmat(1,size(IntraDist)); repmat(2,size(ExtraDist))],1,'+r',0);
	xlabel('Mahalanobis distance')
	title('Box Plot')
	set(gca,'Yticklabel',{'In'; 'Out'},'Xtick',...
		logspace(0,floor(log10(max(Dists))),floor(log10(max(Dists)))+1),'Xlim',[0 max(Dists)],'Xscale','log')
	orient tall
end

if NewDir
	popdir;
end
