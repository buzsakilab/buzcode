function phiG = LinearizeTrajectory(A,ds)
%  
%  USAGE:
%  	phiG = LinearizeTrajectory(A,ds)
%  
%  INPUT:
%  	A: Analysis object
%  	ds: dataset
%  
%  OUTPUT:
%  	phiG: a tsd of linearized positions
%  
%  Adrien Peyrache, 2008



%  The idea here is that phiD/R/L are already normalised. So only the plateform trajectories 
%  have to be normalised, then we "glue" phiD (normalised trajectory on the departure arm) phiP 
%  (the same on plateform) and  phiR/L
%  By putting a negative sign on phiD and by sliding phiD and phiR/L by 0.2 it results in a position which 
%  ranges between -1.2 and 1.2
%  	-1.2 -0.2: departure arm, not always at -1.2 because the animal not systematically goes to the end of 	this arm
%  	-0.2 0.2, plateform
%  	0.2 1.2: arrival arm. Always at least one point at 1.2 because the animal always reaches the end of the arm

A = getResource(A,'PosPhiD',ds);
A = getResource(A,'PosPhiR',ds);
A = getResource(A,'PosPhiL',ds);
A = getResource(A,'PosTimeD',ds);
A = getResource(A,'PosTimeR',ds);
A = getResource(A,'PosTimeL',ds);
A = getResource(A,'PosYS',ds);
A = getResource(A,'PosXS',ds);

phiD = phiD{1};
phiR = phiR{1};
phiL = phiL{1};
XS = XS{1};
YS = YS{1};
timeR = timeR{1};
timeL = timeL{1};
timeD = timeD{1};

rgD = Range(phiD);
rgR = Range(phiR);
rgL = Range(phiL);

% Start time end end time
posSt = min([rgD;rgL;rgR]);
posEn = max([rgD;rgL;rgR]);

% Time on the plateform
timeP = intervalSet(posSt,posEn) - union(timeR,timeL,timeD);
stP = Start(timeP);
enP = End(timeP);

phiP = tsd;

% loop on each trial

for i=1:length(stP)

	% Find first and last coordinates on the plateform for this trial
	xSt = Data(Restrict(XS,stP(i),'align','next'));
	ySt = Data(Restrict(YS,stP(i),'align','next'));
	xEn = Data(Restrict(XS,enP(i),'align','prev'));
	yEn = Data(Restrict(YS,enP(i),'align','prev'));

	% Find extremes points, 
	% x/yD: position just after departure arm,
	% x/yA: position just before arrival arm,  
	
	[yD,ix1] = max([ySt,yEn]);
	[yA,ix2] = min([ySt,yEn]);
	xv = [xSt,xEn];
	xD = xv(ix1);
	xA = xv(ix2);

	% position on the plateform
	x = Restrict(XS,subset(timeP,i));
	y = Restrict(YS,subset(timeP,i));
	pos = [Data(x) Data(y)];

	% v is the vector from entring point to exit point on the plateform
	v = [(xD-xA);(yD-yA)];

	% projects each position position onto v
	p = (pos-repmat([xD,  yD], length(x),1)) * v;

	% normalize p; minimum value:0, max |v|^2
	p = p / norm(v)^2;

	rgP = Range(phiP);
	daP = Data(phiP);

	% concatenation of phiP and the new linearize data
	phiP = tsd([rgP;Range(x)], [daP;p]);

end


rgP = Range(phiP);
pR = Data(phiR);
pL = Data(phiL);
pD = Data(phiD);
pP = Data(phiP);

% normalize pP to range between -0.2 and 0.2
pP = -pP*0.4-0.2;

rg = [rgD;rgL;rgR;rgP];
phi = [-pD-0.2;pL+0.2;pR+0.2;pP];
[rg ix] = sort(rg);
phi = phi(ix);
phiG = tsd(rg,phi);

%  figure(2),clf
%  plot(Range(phiG),Data(phiG))

