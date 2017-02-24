function A = ThetaPfcModulation(A)


A = getResource(A,'PfcThetaPhase');
A = getResource(A,'CellNames');


A = registerResource(A, 'RayleighP_Pfc', 'numeric', {[], []}, ...
    'rayleighP_Pfc', ...
    ['prob of null hypothesis (spike phases are uniformly distr.) based on Rayleigh Z'],'mfile');


A = registerResource(A, 'PrefPh_Pfc', 'numeric', {[], []}, ...
    'prefPh_Pfc', ...
    ['histograms of prefered theta phase of spikes'],'mfile');


A = registerResource(A, 'ModRatio_Pfc', 'numeric', {[], []}, ...
    'modRatio_Pfc', ...
    ['theta modulation ratio of spikes distribution'],'mfile');


resdir = [parent_dir(A), filesep 'PhasePlot'];
[p,ds,e] = fileparts(current_dir(A));

nbCells = length(pfcThetaPhase);

modRatio = zeros(nbCells,1);
rayleighP = zeros(nbCells,1);
prefPh = zeros(nbCells,1);

for j=1:length(pfcThetaPhase)

	phase = Data(pfcThetaPhase{j});
	phase = phase(find(phase>0));
	phase = phase(find(phase<1));
	[modRatio(j), prefPh(j), rayleighP(j), p,H,X,X1]= modRatioFit(2*pi*phase);
	

	fh1 = figure(1),clf
	hold on

	bar(X,H)
	plot(X,p(1)*X1+p(2),'r','LineWidth',2);
	title(['Rayleigh P : ' num2str(rayleighP(j)) ', Mod Ratio : ' num2str(modRatio(j))])

	saveas(fh1,[resdir filesep ds '_' cellnames{j} 'PfcThetaPhase'], 'png');
	

end

if 0
	
	nbPts = size(H,2);
	x = 2*pi/nbPts*[0:nbPts];
	
	figure(1),clf
	hist(modRatio)
	
	figure(2),clf
	hist(prefPh)

end

modRatio_Pfc = {modRatio};
rayleighP_Pfc = {rayleighP};
prefPh_Pfc = {prefPh};

A = saveAllResources(A);

%  keyboard

