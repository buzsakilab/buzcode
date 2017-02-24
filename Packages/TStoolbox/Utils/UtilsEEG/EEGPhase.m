function [ph pk] = EEGPhase(S, eeg, Start, End, minZ)
% ph = EEGPhase(S, EEG, Start, End)
%
% Computes the theta phase of each spike for a group of cells 
%
% INPUTS:
% S:              cell array of ts objects containing the spike trains
% EEG:            EEG tsd object filtered for theta
% Start, end:     vectors containig the start and stop times of the valid
%                 theta epochs
% minZ:		  min z values to detect a peak. EEG has to expressed in zscore.
% OUTPUTS: 
% ph:             cell array of tsd abjects containing the phase of each
%                 spike in S and the theta cycle number

% batta 2000
% peyrache 2007
% status: alpha

%  eeg = Restrict(eeg,intervalSet(Start,End));

dth = diff(Data(eeg));
rg = Range(eeg);

dth1 = [0 dth'];
dth2 = [dth' 0];
clear dth;

t = Range(eeg, 'ts');

peaks = t(find (dth1 > 0 & dth2 < 0));

%  keyboard
clear t;

for iC = 1:length(S)
  s = Range(Restrict(S{iC}, Start, End));
  ph{iC} = zeros(size(s));
  pks = zeros(size(s));
  for j = 1:length(s)
    pk = binsearch_floor(peaks, s(j));

    if pk<length(peaks)
       ph{iC}(j) = (s(j) - peaks(pk)) / (peaks(pk+1) - peaks(pk));
    
       pks(j) = pk;
    end

  end

  ph{iC} = tsd(Range((Restrict(S{iC}, Start, End))), [ph{iC} pks]);
end

pk = Restrict(ts(peaks),intervalSet(s(1),s(end)));
