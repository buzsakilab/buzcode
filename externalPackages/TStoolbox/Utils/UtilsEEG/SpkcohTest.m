Fr1 = 20;
Fr2 = 20;
Tcyc = 1200; % in 10^-4 s, cycle period
TcycE = 50; % in 10^-4 s, error on each period
SpkDisp = 50; %in 10^-4 s, dispersion of spike times
nbCyc = 1000; % 
nbSpkCyc = round(nbCyc/3); % number of cycles where each cells will spikes independantly from one other

T = Tcyc*ones(nbCyc,1)+TcycE*randn(nbCyc,1);

cyc1 = randperm(nbCyc);
cyc1 = sort(cyc1(1:nbSpkCyc));

cyc2 = find(~ismember([1:nbCyc],cyc1));
cyc2 = cyc2(randperm(length(cyc2)));
cyc2 = sort(cyc2(1:nbSpkCyc));

% pour rajouter des spikes dans les memes cycles dans lesquels aucune cellule ne decharge pour l'instant:

%  cyc12 = find(~ismember([1:nbCyc],sort([cyc1 cyc2])));
%  cyc12 = cyc12(randperm(length(cyc12)));
%  cyc12 = cyc12(1:nbSpkCyc/10);
%  cyc1 = sort([cyc1 cyc12]);
%  cyc2 = sort([cyc2 cyc12]);


spk1 = [];

for i=1:length(cyc1)
	nbSpk = poissrnd(Fr1);
	spkTimes = sum(T(1:cyc1(i)))*ones(nbSpk,1)+SpkDisp*randn(nbSpk,1);
	spk1 = [spk1;spkTimes];
end

spk1 = sort(spk1);

spk2 = [];

for i=1:length(cyc2)
	nbSpk = poissrnd(Fr2);
	spkTimes = sum(T(1:cyc2(i)))*ones(nbSpk,1)+SpkDisp*randn(nbSpk,1);
	spk2 = [spk2;spkTimes];
end

spk2 = sort(spk2);

[C,B] = CrossCorr(spk1,spk2,1,500);
figure(1),clf
bar(B,C)

figure(2),clf,hold on
line([spk1';spk1'],(ones(length(spk1),1)*[1 1.3])','Color','r')
line([spk2';spk2'],(ones(length(spk2),1)*[1.3 1.6])','Color','g')

c1 = ts(sort(spk1));
c2 = ts(sort(spk2));
Q = MakeQfromS(tsdArray({c1;c2}),10);
dQ = full(Data(Q));
rg = Range(Q,'s');
Fs = 1./(median(diff(rg)));


movingwin = [10 5];
params.fpass = [0 20];
params.tapers = [3 5];
params.Fs = Fs;
params.err = [1 0.05];
[C,phi,S12,S1,S2,t,f]=cohgramc(dQ(:,1),dQ(:,2),movingwin,params);

figure(3),clf
Csm = convn(C,gausswin(6),'same');
imagesc(t(3:end-3),f,Csm'),axis xy

C12 = corrcoef([dQ(:,1) dQ(:,2)])
