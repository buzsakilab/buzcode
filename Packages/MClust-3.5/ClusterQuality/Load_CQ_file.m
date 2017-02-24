function stats = Load_CQ_file(fc)

% stats = Load_CQ_file(fc);
% 
% For a cell array of strings listing .t files to get, loads the
% appropriate -ClusterQual.mat file and puts the values for each .t file
% into an output structure, which has fields: nSpikes, L (amount of
% contamination), L_ratio (L/nSpikes), IsolationDist, and SNR
% (signal-to-noise ratio, defined by GetSNR).

% ncst 02 Dec 04

nFiles = length(fc);

stats.fc = fc;
stats.nSpikes = repmat(nan,size(fc));
stats.L = repmat(nan,size(fc));
stats.L_Ratio = repmat(nan,size(fc));
stats.IsolationDist = repmat(nan,size(fc));
stats.SNR = repmat(nan,size(fc));

for iF = 1:nFiles
	Name = FindNameInfo(fc(iF));
	CQ_target_fn = [Name.Location Name.CellID '-ClusterQual.mat'];
	if exist(CQ_target_fn)
		
		temp = load(CQ_target_fn);
		Fields_temp = fields(temp);
		
		if ~isempty(strmatch('CluSep',Fields_temp))
			temp = temp.CluSep;
			Fields_temp = fields(temp);
			
			% check to make sure the appropriate fields exist, and if so, then
			% load the value into the output structure
			if ~isempty(strmatch('nSpikes',Fields_temp))
				stats.nSpikes(iF) = temp.nSpikes;
			end
			if ~isempty(strmatch('L',Fields_temp))
				stats.L(iF) = temp.L;
			end
			if ~isempty(strmatch('Lratio',Fields_temp))
				stats.L_Ratio(iF) = temp.Lratio;
			end
			if ~isempty(strmatch('IsolationDist',Fields_temp))
				stats.IsolationDist(iF) = temp.IsolationDist;
			end
			if ~isempty(strmatch('SNR',Fields_temp))
				stats.SNR(iF) = max(temp.SNR);
			end
		end
	else
		disp(['  ' Name.CellID '-ClusterQual.mat file not found'])
	end
	
	displayprogress(iF,nFiles)
end