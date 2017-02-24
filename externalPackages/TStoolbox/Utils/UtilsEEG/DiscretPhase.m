function ph = DiscretPhase(S, tsa, ep)
% ph = DiscretPhase(S, CRtheta, ThStart, ThEnd)
%
% Computes the theta phase of each spike for a group of cells 
%
% INPUTS:
% S:              cell array of ts objects containing the spike trains
% tsa:        	  a tsd of events the phase is calculated with. each event represents the 0 phase of a underlying oscillation.
% OPTION:
% ep:		  a intervalSet defining the time the phase has to ba computed.
% OUTPUTS: 
% ph:             cell array of tsd abjects containing the phase of each
%                 spike in S and the cycle number

% batta 2000
% peyrache 2007
% status: alpha

peaks = Range(tsa);

if nargin<3

	ep = timeSpan(tsa);

end

for iC = 1:length(S)
	s = Range(Restrict(S{iC}, ep));
	ph{iC} = zeros(size(s));
	pks = zeros(size(s));

	for j = 1:length(s)
		pk = binsearch_floor(peaks, s(j));
		ph{iC}(j) = (s(j) - peaks(pk)) / (peaks(pk+1) - peaks(pk));
		pks(j) = pk;
	
	end
	
	ph{iC} = tsd(s, [ph{iC} pks]);
end
