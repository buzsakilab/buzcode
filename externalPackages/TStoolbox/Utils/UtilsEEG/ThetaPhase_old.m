function ph = ThetaPhase(S, CRtheta, ThStart, ThEnd)
% ph = ThetaPhase(S, CRtheta, ThStart, ThEnd)
%
% Computes the theta phase of each spike for a group of cells 
%
% INPUTS:
% S:              cell array of ts objects containing the spike trains
% CRtheta:        EEG tsd object filtered for theta
% ThStart, Thend: vectors containig the start and stop times of the valid
%                 theta epochs
%
% OUTPUTS: 
% ph:             cell array of tsd abjects containing the phase of each
%                 spike in S and the theta cycle number

% batta 2000
% peyrache 2007
% status: alpha



dth = diff(Data(CRtheta));

dth1 = [0 dth'];
dth2 = [dth' 0];
clear dth;

t = Range(CRtheta, 'ts');

thpeaks = t(find (dth1 > 0 & dth2 < 0));
clear t;
ph = cell(length(S),1);
for iC = 1:length(S)
  s = Range(Restrict(S{iC}, ThStart, ThEnd));
  ph{iC} = zeros(size(s));
  pks = zeros(size(s));
  for j = 1:length(s)
    pk = binsearch_floor(thpeaks, s(j));
	
    t = (s(j) - thpeaks(pk)) / (thpeaks(pk+1) - thpeaks(pk));
    ph{iC}(j) = 2*pi*mod(t,1);
    pks(j) = pk;
    
  end

%    ph{iC} = tsd(Range((Restrict(S{iC}, ThStart, ThEnd))), [ph{iC} pks]); %);
  ph{iC} = tsd(Range((Restrict(S{iC}, ThStart, ThEnd))), ph{iC}); %);

end