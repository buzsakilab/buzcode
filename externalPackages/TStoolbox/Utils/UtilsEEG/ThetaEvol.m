dset = 'Rat18/181011';

A = getResource(A,'RatioRipS2_2',dset);
A = getResource(A,'RatioRipS1_2',dset);
A = getResource(A,'TrialOutcome',dset);
A = getResource(A,'ProbaTheta',dset);

A = getResource(A,'HcThetaPhase',dset);
A = getResource(A,'HcTrace',dset);
st = Range(trialOutcome{1});
hcThetaPhase = hcThetaPhase{hcTrace-4};
trials = intervalSet(st(1:end-1),st(2:end));
modThetaMat=[];
prefPhMat = [];
for i=1:length(hcThetaPhase)
for j=1:length(st)-1
ph = 2*pi*Data(Restrict(hcThetaPhase{i},subset(trials,j)));
if length(ph)>10
[prefPhMat(i,j), Rmean, delta, p] = CircularMean(ph);
[modThetaMat(i,j) b c d]= modRatioFit(ph,prefPhMat(i,j));
end;
end
end

[m ix] = sort(probaTheta);
imagesc([log(ratioRipS1_2(ix)) (modThetaMat(ix,:)) log(ratioRipS2_2(ix))]);