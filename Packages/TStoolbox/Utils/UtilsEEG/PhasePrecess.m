function A = PhasePrecess(A)

%  A = getResource(A,'PfcThetaPhase');
A = getResource(A,'HcThetaPhase');

A = getResource(A,'StartTrial')
startTrial = startTrial{1};

A = getResource(A,'TrialOutcome')
trialOutcome = trialOutcome{1};

%  A = getResource(A,'PosYS');
%  YS = YS{1};

%  A = getResource(A,'PosPhiR');
%  phiR = phiR{1};
%  
%  A = getResource(A,'PosPhiL');
%  phiL = phiL{1};
%  
%  A = getResource(A,'PosPhiD');
%  phiD = phiD{1};
%  
%  a = getResource(A,'PosOffRewardL');
%  offRewardL = offRewardL{1};
%  
%  A = getResource(A,'MazeEpoch');
%  mazeEpoch = mazeEpoch{1};

%  A = getResource(A,'SpikeData');
%  
%  Q = MakeQfromS(S,1ss500);
%  Q = Restrict(Q,mazeEpoch{1});
%  dQ = Data(Q);
%  Qr = tsd(Range(Q),dQ(:,32));

st = Range(startTrial)-10000;
to = Range(trialOutcome);
side = Data(trialOutcome);

%  offRewardL = threshold(Data(phiL), 0.90, 'Crossing', 'Falling');

st = st(side==0);
to = to(side==0);

%  trials = mazeEpoch;
trials = intervalSet(st,to);
%  trials = intervalSet(Restrict(phiL,st,'align','next'),to);
%  trials = intervalSet(to,Restrict(offRewardL,to,'align','next'));
%  trials = intervalSet(Restrict(offRewardL,to,'align','next'),Restrict(phiD,to,'align','next'));


%  pfcPh = Restrict(pfcThetaPhase{1},trials);
%  hcPh = Restrict(hcThetaPhase{1},trials);

%  
%  figure(1),clf
%  hold on
%  plot(Data(Restrict(phiL,pfcPh)),Data(pfcPh),'ko');
%  plot(Data(Restrict(phiL,pfcPh)),Data(pfcPh)-1,'ko');
%  
%  
%  figure(2),clf
%  hold on
%  plot(Data(Restrict(phiL,hcPh)),Data(hcPh),'ko');
%  plot(Data(Restrict(phiL,hcPh)),Data(hcPh)-1,'ko');

%  
%  
%  A = getResource(A,'PosYS');
%  A = getResource(A,'PosXS');
%  
%  keyboard
%  
%  YS = Restrict(YS{1},trials);
%  XS = Restrict(XS{1},trials);
%  
%  
%  
%  %  keyboard 
%  PhaseDerivMap(XS,YS,hcPh);
%  PhaseMap(XS,YS,hcPh);
%  
%  keyboard


figure(3),clf
hold on;

for i=1:length(to)

%  	pfcPh = Restrict(pfcThetaPhase{1},offRrew(i),offRarm(i));
	
	hcPh = Restrict(hcThetaPhase{48},subset(trials,i));
	l = to(i)-st(i);
	plot((Range(hcPh)-st(i))/l,Data(hcPh),'k.')
	plot((Range(hcPh)-st(i))/l,Data(hcPh)-1,'k.')	

end

%  
%  figure(1),clf
%  plot(Data(Restrict(Qr,pfcPh)),Data(pfcPh),'ko');
%  
%  figure(2),clf
%  plot(Data(Restrict(Qr,hcPh)),Data(hcPh),'ko');