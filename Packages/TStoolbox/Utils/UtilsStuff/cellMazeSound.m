A = getResource(A,'PosYS','Rat18/181014');
A = getResource(A,'PosXS','Rat18/181014');
XS = XS{1};
YS = YS{1};

A = getResource(A,'SpikeData','Rat18/181014');
A = getResource(A,'MazeEpoch','Rat18/181014');
mazeEpoch = mazeEpoch{1};


tSpk = Range(Restrict(S{32},mazeEpoch));

XS = Restrict(XS,mazeEpoch);
YS = Restrict(YS,mazeEpoch);

pX = Restrict(XS,tSpk);
pY = Restrict(YS,tSpk);

tBegin = Start(mazeEpoch);
maxX = max(Data(XS));
minX = min(Data(XS));
maxY = max(Data(YS));
minY = min(Data(YS));


stTime = Start(mazeEpoch)+30000;

figure(1),clf

for i=stTime+900000:300:stTime+2000000%End(mazeEpoch)
	
	clf
%  	display(num2str(i-stTime))
	hold on
	plot([minX minX maxX maxX],[minY maxY maxY minY],'k*')
	plot(Data(Restrict(XS,i-30000,i)), Data(Restrict(YS,i-30000,i)))
	plot(Data(Restrict(pX,i-30000,i)), Data(Restrict(pY,i-30000,i)),'ro')
	
	hold off
	

	pause(0.01);

end



