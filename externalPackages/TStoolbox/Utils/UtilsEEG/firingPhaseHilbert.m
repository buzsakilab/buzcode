function ph = firingPhaseHilbert(eeg, S,W,epoch) 
% ph = thetaPhaseHilbert(CRtheta,S) copmute firing phase of spikes with the Hilbert Transform 
%
% INPUTS:
% CR: tsd containing EEG filtered in the frequency domain of interest (e.g.
%        theta, or gamma
% S (optional): tsd array of cells, the firing phase of each spike will be calculated 
% W: pass-ban vector (normalized between 0 and 1 to be used by cheby2)
% OUTPUT: 
% ph: a tsdArray of firing phases for the spike trains in S

% copyright (c) 2006-2011 Francesco P. Battaglia & Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
%version 0.1 under construction

[b,a] = cheby2(4,40,W);

dEeg = Data(eeg);
dEeg = filtfilt(b,a,dEeg);
dEeg = hilbert(dEeg);
dEeg = mod(atan2(imag(dEeg), real(dEeg)),2*pi);
eeg = tsd(Range(eeg),dEeg);

eeg = Restrict(eeg,epoch);
S = Restrict(S,epoch);

ph = cell(length(S),1);


for ii = 1:length(S)
	ph{ii} = Restrict(eeg, S{ii}, 'align', 'closest', 'time', 'align');
end

ph = tsdArray(ph);