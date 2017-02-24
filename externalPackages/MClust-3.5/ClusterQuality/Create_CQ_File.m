function Create_CQ_File(fc,recalc,fd_TT)

% Create_CQ_file(fc,recalc,fd_TT);
% 
% Used to create a -ClusterQual.mat file for .t spike train files
% 
% Inputs: fc - cell array of strings indicating each .t file to create.
%              Assumes .t files are named with standard Rat SSN conventions
%              (i.e. RXXX-YYYY-MM-DD-TTXX-XX.t).  If fc is empty, will find
%              all .t files in the current directory.
%         recalc - if 1, even if a cluster quality file exists, the values
%                  will be recalculated, default = 0.
%         fd_TT - function assumes the tetrode files for each .t file are
%                 located in the same directory as the .t file.  If not
%                 (like when files are on CD) you can specify a directory
%                 to check, like a CD drive
%
% Outputs: Saves a -ClusterQual.mat file, with the stem being the Rat SSN
%          and tetrode/cluster designations.  File contains a structure,
%          CluSep, the fields of which are the IsolationDistance, L,
%          L_ratio, number of spikes, features used to calculate quality,
%          and the degrees of freedom (dimensionality of featurespace).
%          Signal-to-noise ratio is also saved for each waveform channel.

% ncst 02 Dec 04
% modified ADR Feb 2008

% replace with formats you use. _ID_ will be replaced with the tetrode number
FileFormats = {'R*TT_ID_.ntt'};
%FileFormats = {'R*TT_ID_.ntt', 'TT_ID_.ntt', 'R*TT_ID_.dat*','R*Sc_ID_.ntt*', 'TT_ID_.dat*' 'Sc_ID_.ntt*', 'TT_ID_.tt*', 'TT0_ID_.tt*'};


% initialize variables which were not passed in.
if ~exist('fc','var'); fc = FindFiles('*.t'); end;
if ~exist('recalc','var'); recalc = 0; end;
if ~exist('fd_TT','var'); fd_TT = ' '; end;

% Assume this is hard coded (the feature space to use is always the same)
FeaturesToUse = {'energy','wavePC1'};
nRecords_SNR = 20000; % maximum number of records to load when calculating the signal-to-noise ratio
				
if isempty(fc); fc = FindFiles('*.t'); end;

fc_remove = {}; % any tetrode files that are copied from CD and need to be erased.
fd_Curr = ' '; % current directory (of last tetrode loaded)
TT_Curr = ' '; % current tetrode (last one loaded)

nAttempted = 0; % number of .t files attempted
nGood = 0;     % number of .t files successfully processed

nFeatures = length(FeaturesToUse);
S = LoadSpikes(fc,'tsflag','ts'); % load the spikes for each .t file

for iFC = 1:length(fc)
	recalcYN = recalc;  % start out with whatever our passed in value was.
    Name = FindNameInfo(fc{iFC});
	
	fc_CQ = [Name.CellID '-ClusterQual.mat'];  % file to be created

	disp(' ');
	disp(['Processing file ' num2str(iFC) ' of ' num2str(length(fc)) ' ' Name.CellID])
    pushdir(Name.Location);
	
	% Check the existing ClusterQuality.mat file to see if we need to
	% recalculate it
	clear CluSep
	if exist(fc_CQ,'file') && recalcYN == 0 % if there is a CQ file
		load(fc_CQ)
		if exist('CluSep','var') % Check to see if it is recent enough to use the CluSep structure
			if ~isempty(strmatch('Features',fields(CluSep))) % Check for a features field
				if length(CluSep.Features) ~= length(FeaturesToUse) % Check for the correct number of features
					recalcYN = 1;
				end
				
				for iF = 1:length(FeaturesToUse) % Check each feature to see if it is one we want to use
					if isempty(strmatch(FeaturesToUse{iF},CluSep.Features))
						recalcYN = 1;
					end
				end
			end
		else
			recalcYN = 1;
		end
	else
		recalcYN = 1;
	end
	
	if ~recalcYN % if we should try to calculate the CQ file	
		disp('  Using existing file...')
	else 
		nAttempted = nAttempted + 1; 
		removeYN = 0;
		Tetrode = num2str(str2num(Name.Tetrode), '%02d');  %#ok<ST2NM>
	
		% Look for the tetrode using known tetrode naming formats
		fc_TT = {};
		for iFT = 1:length(FileFormats)
			fc_TT = cat(1, Findfiles(strrep(FileFormats{iFT}, '_ID_', Tetrode)));
		end
		
		% if we did not find a TT file
		if isempty(fc_TT) && isempty(strmatch(fd_TT,' ')) % if a directory was supplied to check for the tetrode file (like from a CD) move to that dir
			pushdir(fd_TT);
			for iFT = 1:length(FileFormats)
				fc_TT = cat(1, Findfiles(strrep(FileFormats{iFT}, '_ID_', Tetrode)));
			end
			popdir;
	
% 			if ~isempty(fc_TT)
% 				dos(['copy ' fc_TT{1} ' ' Name.Location]);
% 				removeYN = 1;
% 			end
		
		end
	
		if isempty(fc_TT); % if we still do not have a tetrode file, we cant' create the CQ file
			disp(['No tetrode file found for ' Name.CellID ]);
		else
			%-------------------------------------------------------------
			% Check to see if the tetrode file needs to be unzipped
			[p n e] = fileparts(fc_TT{1});
			isZipped = strcmpi('.gz',e);
			if isZipped
				unix(['gunzip ',  char(fc_TT(1))]);
				
				fc_TT = [findfiles(['TT' Tetrode '.dat*']) findfiles(['Sc' Tetrode '.ntt*']) ...
						findfiles(['TT' Tetrode '.tt*']) findfiles(['TT0' Tetrode '.tt*'])];
			end
			
			%-------------------------------------------------------------
			if removeYN; fc_remove(end+1) = fc_TT(1); end; % if we pulled this file off of a CD
				
			% get the mclust feature data file defaults, set the maximum
			% number of records for .ntt and .dat files (which are smaller
			% than the 64 sample waveforms I've used from 32kHz recordings
			Get_MClust_FD_Defaults;
			[p n e] = fileparts(fc_TT{1});
			if strcmpi('.ntt',e) || strcmpi('.dat',e)
				record_block_size = 100000;
			end
			
			if ~strcmpi(Name.Location,fd_Curr) || ~strcmpi(Name.Tetrode,TT_Curr) % if the tetrode to load is not the one currently in memory
				nS = MClust_CountSpikes(fc_TT{1});  % find the total number of records
				if nS < record_block_size  % if the number of records is smaller than our maximum, get all of the records
					[T WVD] = MClust_LoadNeuralData(fc_TT{1});
				else % if not, get all of the timestamps
					T = MClust_LoadNeuralData(fc_TT{1});
					WVD = [];
				end
			end
		
			badTT = 0;  % Are the spikes aligned incorrectly?
			ds = Data(S{iFC});
			f_cluster = interp1([(T(1) - max(diff(T))); T; T(end) + max(diff(T))],0:(nS + 1),ds,'nearest');
			f_cluster(f_cluster > nS) = NaN;
			if any(isnan(f_cluster))
				disp(['Lost spikes for ' Name.SessionID])
				beep
				badTT = 1;
			end
			
			% check to see if the records we aligned to are within 1
			% timestamp of the times in the .t file, if not, we have loaded
			% the wrong tetrode OR the wrong loading engine
			if ~badTT 
				if max(abs(T(f_cluster) - ds)) > 1
					disp(['Timestamp mismatch for ' Name.SessionID])
					beep
					badTT = 1;
				end
			end
	
			if badTT % if we got the wrong data loaded
				disp(['*** Bad or missing tetrode file for ' Name.CellID]);
			else
				if length(f_cluster) > nRecords_SNR
					X = randperm(length(f_cluster));
					f = f_cluster(X(1:nRecords_SNR));
				else
					f = f_cluster;
				end
				
				if isempty(WVD)
					[t,WV] = MClust_LoadNeuralData(fc_TT{1},f,2);
				else
					t = T(f);
					WV = WVD(f,:,:);
				end
				
				[jnk nCh nWVSamples] = size(WV);
				
				%-------------------------------------------------------------
				% While we are at it, check for the existence of the average
				% waveform file, if it doesn't exist, create it
				fc_WV = [Name.Location Name.SessionID '-TT' Name.Tetrode '-' Name.Cluster '-wv.mat'];
				if ~exist(fc_WV,'file')
					disp('  Creating average waveform file')
					for it = 1:4
						mWV(:,it) = squeeze(mean(WV(:,it,:),1));
						sWV(:,it) = squeeze(std(WV(:,it,:),1));
						xrange(:,it) = (((nWVSamples + 2) * (it-1)) + (1:nWVSamples))';
					end
					save(fc_WV,'mWV','sWV','xrange','-mat');
				end
				
				%-------------------------------------------------------------
				% Signal to noise ratio calculation
				disp('  Calculating signal-to-noise ratio')
				ClustTT = tsd(t,WV);
				
				if nS > nRecords_SNR
					X = randperm(nS);
					f = X(1:nRecords_SNR);
				else
					f = 1:nS;
				end
				
				if isempty(WVD)
					[t,WV] = MClust_LoadNeuralData(fc_TT{1},f,2);
				else
					t = T(f);
					WV = WVD(f,:,:);
				end
				NoiseTT = tsd(T,WV);
				
				SNR = GetSNR(1,'ClustTT',ClustTT,'NoiseTT',NoiseTT);
				
				%-------------------------------------------------------------
				% Cluster quality
				disp('  Calculating cluster quality')
				
				% Create feature space
				if ~strcmpi(Name.Location,fd_Curr) || ~strcmpi(Name.Tetrode,TT_Curr)
					if isempty(WVD)
						Fet = repmat(NaN,nS,nCh*nFeatures);
						nB = ceil(nS/record_block_size);
						for iB = 1:nB
							f_get = (1:record_block_size) + (iB - 1)*record_block_size;
							f_get = f_get(f_get <= nS);
							
							[t,WVD] = MClust_LoadNeuralData(fc_TT{1},f_get,2);
							
							Fet(f_get,:) = Create_FeatureSpace(WVD);
						end
						WVD = [];
					else
						Fet = Create_FeatureSpace(WVD);
					end
					
					fd_Curr = Name.Location;
					TT_Curr = Name.Tetrode;
				end
				
				CluSep = Cluster_Quality(Fet,f_cluster);
				CluSep.SNR = SNR;
				CluSep.Features = FeaturesToUse;
				CluSep.nSpikes = length(ds);
				CluSep.nEvents = nS;
				
				save(fc_CQ,'CluSep'); % Save output file
				nGood = nGood + 1;
			end
		end
	end
end
		
for iFC = 1:length(fc_remove)
	unix(['del ',  char(fc_remove(iFC))]);
end

disp(' ');
disp([' ' num2str(nGood) '/' num2str(nAttempted) ' files completed successfully.']);
disp(' ');