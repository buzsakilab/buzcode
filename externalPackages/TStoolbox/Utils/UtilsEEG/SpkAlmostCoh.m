Fr1 = 20;
Fr2 = 20;
Tcyc = 1200; % in 10^-4 s, cycle period
TcycE = 0; % in 10^-4 s, error on each period
SpkDisp = 0; %in 10^-4 s, dispersion of spike times
nbCyc = 1200; % 
nbSpkCyc = round(nbCyc/3); % number of cycles where each cells will spikes independantly from one other
epsilon = 1;

T = Tcyc*ones(nbCyc,1)+TcycE*randn(nbCyc,1);
t = [0:nbCyc*Tcyc];
s = sin(2*pi*t/Tcyc);
cyc1 = randperm(nbCyc);
cyc1 = sort(cyc1(1:nbSpkCyc));


spk1 = cumsum(T-epsilon);
phase = (spk1-cumsum(T))./Tcyc;
[C,B] = CrossCorr(spk1,spk1,1,500);
figure(1),clf
bar(B,C)

figure(2),clf,hold on
plot(t,s)
line([spk1';spk1'],(ones(length(spk1),1)*[1 1.3])','Color','r')

figure(3),clf
hist(phase,100)

c1 = ts(sort(spk1));
Q = MakeQfromS(tsdArray({c1}),1);
dQ = full(Data(Q));
rg = Range(Q,'s');
Fs = 1./(median(diff(rg)));


movingwin = [10 5];
params.fpass = [0 20];
params.tapers = [3 5];
params.Fs = Fs;
params.err = [1 0.05];
[C,phi,S12,S1,S2,t,f]=cohgramc(s(1:1437602)',dQ(:,1),movingwin,params);

figure(3),clf
%  Csm = convn(C,gausswin(6),'same');
imagesc(t(3:end-3),f,C'),axis xy

C12 = corrcoef([dQ(:,1) dQ(:,2)])
