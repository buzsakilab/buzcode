function [cmats_sig,strength,pref_phase,mean_phase] = thetamod2(LFPtheta,signals,freqs,phaseBand,Fs,band)
%   Detailed explanation goes here

%% initializa
t = (0:(size(LFPtheta,1)-1))'/Fs;
% Find theta epochs  

% Theta parameters using Hilbert transform
theta = filtsig(LFPtheta, 1000/Fs, [phaseBand(1) phaseBand(2) phaseBand(3) phaseBand(4)]/1000);
%theta = filtsig(LFPtheta, 1000/Fs, [2 4 12 14]/1000); default for theta
hilb = hilbert(theta);
thetaP = angle(hilb);
%thetaH = abs(hilb);

% Theta parameters using Belluscio's method
%[thetaP, thetaH] = thetaParams_fromExtrema (LFP, Fs, thetaEpochInds);

% Theta phase bins
nbins_theta = 32;
xtheta = 0:(2*pi/nbins_theta):(2*pi);
[~,thetabins] = histc(mod(thetaP, 2*pi), xtheta);

LFP_thetaRUNAvg = zeros(nbins_theta,size(LFPtheta,2));
for i=1:nbins_theta
     LFP_thetaRUNAvg(i,:) = mean(LFPtheta((thetabins == i),:));
end

%% Spectrograms
fspec = freqs';
cwtScales = Fs./fspec;

signalsCWTs = zeros(length(fspec), size(signals,1), size(signals,2));
for i=1:size(signals,2)
    tmpcwt = abs(cwt(signals(:,i), cwtScales, 'cmor3-1'));
    signalsCWTs(:,:,i) = tmpcwt;
end
clear tmpcwt
ThetaSpectra = zeros(length(fspec), nbins_theta, size(signals,2));

for i=1:nbins_theta
    ThetaSpectra(:,i,:) = mean(signalsCWTs(:,(thetabins == i),:), 2);
end

%% Plot theta modulation
plotsignals = 1:size(signals,2);
plotfreqs = 1:length(fspec);

% Heat map matrices
cmats_sig = zeros(length(plotfreqs), nbins_theta, length(plotsignals),2);
for j=1:length(plotsignals)
    for k=1:length(plotfreqs)
        cmats_sig(k,:,j,1) = zScore(ThetaSpectra(plotfreqs(k),:,j), mean(signalsCWTs(plotfreqs(k),:,j)), ...        
             std(signalsCWTs(plotfreqs(k),:,j)));
    end
end
clims = max(maxall(abs(cmats_sig)))*[-1 1];

figure
axs = zeros(length(plotsignals),2);
imgs = zeros(length(plotsignals),2);
lines = zeros(length(plotsignals),2);
axs(1,1) = subplot(1,length(plotsignals),1);
for j=1:length(plotsignals)
    axs(j,1) = subplot(1,length(plotsignals),j);
    imgs(j,1) = imagesc([midpoints(xtheta(1:2)) 4*pi-midpoints(xtheta(1:2))]*180/pi, fspec(plotfreqs([1 end])), cmats_sig(:,[1:end 1:end],j,1), clims);
    set(gca, 'YDir', 'normal', 'XLim', [0 720], 'XTick', 0:180:720);
    hold on;
    lines(j+1,1) = plot([midpoints(xtheta) midpoints(xtheta)+2*pi]*180/pi, (cos([midpoints(xtheta) midpoints(xtheta)]) - 1)*10 + 150, 'k--');
    hold off;
if j>1
    set(gca,'YTick',[]);
end
end

%% cuantifications
phases=0:11.612:360; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% poner esto en función de thetabins no asi 
fc1v=abs(freqs-band(1));
fc1= find (fc1v==min(fc1v)); 
fc2v=abs(freqs-band(2));
fc2= find (fc2v==min(fc2v)); 
tetmod=cmats_sig(fc1:fc2,:,:,:);
% modulation stregth
strength=zeros(size(signals,2),1);
for i=1:size(signals,2)
    for j=1:length(tetmod(:,1,1,1))
        aux(j,i)=max(tetmod(j,:,i,1));
    end
    strength(i,1)=mean(aux(:,i));
end
% prefered phase
pref_phase=zeros(size(signals,2),1);
aux3=zeros(length(tetmod(:,1,1,1)),size(signals,2));
for i=1:size(signals,2); 
    for j=1:length(tetmod(:,1,1,1)) % n frecs selected
        aux(j,i)=max(tetmod(j,:,i,1));
        aux2(j,i)=find(tetmod(j,:,i,1)==aux(j,i));
        aux3(j,i)=phases(aux2(j,i));
    end
        pref_phase(i,1)= mean (aux3(:,i));
end
% Mean phase 
pow_phase=zeros(size(tetmod(1,:,1,1),size(signals,2))); 
for i=1:size(signals,2)
    for j=1:length(tetmod(1,:,1,1)) % para cada bin de fase
        pow_phase(j,i)=mean(tetmod(:,j,i,1)); % pow_phase is the phase distrib of power
    end
end 
for i=1:size(pow_phase,2)
    z=pow_phase(:,i)'.*exp(1i*phases);
    n=sum(pow_phase(:,i));
    zm=sum(z)./n;
    mean_phase(i)=abs(angle(zm));
end
mean_phase=(mean_phase'.*180)/pi;
end