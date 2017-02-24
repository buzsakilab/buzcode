function A = RosePlot(A);

A = getResource(A,'HcThetaPhase');
A = getResource(A,'CellNames');


for j=32:32% length(hcThetaPhase)

	phase = Data(hcThetaPhase{j});
	phase = phase(find(phase>0));
	phase = phase(find(phase<1));
	[modRatio(j), prefPh(j), rayleighP(j), p,H,X,X1]= modRatioFit(2*pi*phase);
	

	fh1 = figure(1),clf
	hold on

	bar(X,H)
	plot(X,p(1)*X1+p(2),'r','LineWidth',2);
	title(['Rayleigh P : ' num2str(rayleighP(j)) ', Mod Ratio : ' num2str(modRatio(j))])

	saveas(fh1,[resdir filesep ds '_' cellnames{j} 'ThetaPhase'], 'png');
	

	phase = Data(pfcThetaPhase{j});
	phase = phase(find(phase>0));
	phase = phase(find(phase<1));
	[modRatio(j), prefPh(j), rayleighP(j), p,H,X,X1]= modRatioFit(2*pi*phase);
	

	fh1 = figure(1),clf
	hold on

	bar(X,H)
	plot(X,p(1)*X1+p(2),'r','LineWidth',2);
	title(['Rayleigh P : ' num2str(rayleighP(j)) ', Mod Ratio : ' num2str(modRatio(j))])

	saveas(fh1,[resdir filesep ds '_' cellnames{j} 'ThetaPhase'], 'png');
	


end
