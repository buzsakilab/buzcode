function [modRatio,p,H,X,X1]= modRatioFit(phase,prefPh)

%  [modRatio, prefPh, fitCurve, H]= modRatioFit(phase) 
%   Adrien Peyrache 2007


nbSpk = length(phase);
[H,X] = hist(phase,2*pi*[0:0.01:1]);
X1 = cos(X-prefPh);
[p S] = polyfit(X1,H,1);
modRatio = abs(p(1)/p(2));


if 0

	figure(3),clf
	hold on

	bar(X,H)
	plot(X,p(1)*X1+p(2));

end
%  keyboard