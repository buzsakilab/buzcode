function KKStatistics(iC)

% KKStatistics(iC)
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.


   global MClust_TTData MClust_FeatureSources
   global MClust_FeatureIndex %ncst
   global KKClust
   global KlustaKwik_Clusters
   
% Statistics text
f = FindInCluster(KlustaKwik_Clusters{iC});
[TT,was_scaled] = ExtractCluster([], f);
msgstr = Stats(iC, TT,was_scaled);
KKClust.Stats{iC} = msgstr;
  
% Average Waveforms
if length(Range(TT,'ts')) > 1
    [wfm wferr] = AverageWaveform(TT);
    KKClust.WaveForms{iC} = {wfm, wferr};
else
    wfm = squeeze(Data(TT));
    KKClust.WaveForms{iC} = {wfm, zeros(size(wfm))};
end
 
% ISI histograms
[H binsUsed] = HistISI(ts(Range(TT,'ts')));
KKClust.ISI{iC} = {H, binsUsed};


%=============================================================
function msgstr = Stats(iC, TT,was_scaled)
global KlustaKwik_Clusters % added number of clusters to message string -- JCJ Sept 2007

msgstr = {['Cluster ' num2str(iC) ' of ' num2str(length(KlustaKwik_Clusters))]}; % added number of clusters to message string -- JCJ Sept 2007
if was_scaled == 0
    nSpikes = length(Range(TT, 'ts'));
else
    nSpikes = was_scaled;
end

if nSpikes > 1
	msgstr{end+1} = sprintf('%.0f spikes ', nSpikes);
	dts = diff(Range(TT, 'ts'));
	
	mISI = mean(dts)/10000;
	msgstr{end+1} = sprintf('inv mean ISI = %.4f spikes/sec ', 1/mISI);
	mdISI = median(dts)/10000;
	if isempty(mdISI)
        mdISI = 0;
	end
	mdFR = 1/mdISI;
	msgstr{end+1} = sprintf('inv med. ISI = %.4f spikes/sec ', mdFR);
	fr = 10000 * nSpikes/(EndTime(TT) - StartTime(TT));
	msgstr{end+1} = sprintf('firing rate = %.4f spikes/sec ', fr);
else
    msgstr{end + 1} = 'Only 1 spike in this cluster';
end