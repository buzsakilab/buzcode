function PhaseMap = PhaseMap(X,Y,ph)

%peyrache 2007, beta

times = Range(ph);
phase = Data(ph);

[occH, xbin, ybin] = hist2d(Data(X), Data(Y), 200, 200);

PhaseMap = zeros(length(xbin),length(ybin));
phH = zeros(length(xbin),length(ybin));

for j=1:length(ph)

	t = times(j);
	x = Data(Restrict(X,t));
	y = Data(Restrict(Y,t));
	xi = binsearch(xbin,x);
	yi = binsearch(ybin,y);

	PhaseMap(xi,yi) = PhaseMap(xi,yi) + phase(j);
	phH(xi,yi) = phH(xi,yi) + 1;
%  	display([num2str(xi) ',' num2str(yi) ': ' num2str(phase(j))])

end  

dbclear warning
warning off
PhaseMap = PhaseMap./phH;
warning on
PhaseMap(isnan(PhaseMap))=0;
largerMatrix = zeros(250,250);
largerMatrix(26:225,26:225) = PhaseMap;
PhaseMap = largerMatrix;
warning on


for x=1:19;
	gw(x) = exp(-((x-10)*(x-10))/(40));
end

gw = diff(gw);

PhaseMap = convn(PhaseMap,gw,'same')/sum(gw);
PhaseMap = convn(PhaseMap',gw,'same')/sum(gw);


%  for x=1:19;for y=1:19; gw(x,y) = exp(-((x-10)*(x-10) + (y-10)*(y-10))/(2*20));end;end
%  PhaseMap = conv2(PhaseMap,gw,'same')/sum(sum(gw));
%  
%  figure(2),clf
%  imagesc(phH)
%  colorbar

figure(1),clf
imagesc(PhaseMap)
colorbar

%  keyboard