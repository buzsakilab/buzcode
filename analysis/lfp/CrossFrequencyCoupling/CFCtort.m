function [tort,freq,MI] = CFCtort(LFPtheta,signals,frecs1,frecs2,Fs,band)
% [tort,freq,MI] = CFCtort(LFPtheta,signals,freqs1,freqs2,Fs,band)
% Computes and plots phase - amplitude cross-frequency coupling and calculates de
% modualtion index (Tort et al., 2009)
%
% Inputs:   LFPtheta = LFP signal to extract the phase (only one)
%           signals  = LFP signals to extract the amplitude
%           freqs1   = vector with the frequency values to extract the
%                      phase (e.g.: freqs1= 2:0.25:15;)
%           freqs2   = vector with the frequency values to extract the
%                      amplitude (e.g.: freqs1= 30:2:150;)
%           Fs       = sampling rate in Hz
%           band     = frequency band of interest to calculate the MI
% Outputs:  tort     = phase - amplitude comodulogram
%           freq     = frequency of maximun amplitude modulation (preferred frequency)
%           MI       = modulation index

%%%% Antonio Fern�ndez-Ruiz, 2014.

% THIS IS A TEST OF GIT FUNCTIONALITY

%% 
Nptos=size(signals,1);
Nsenalesamp=size(signals,2); 
Nsenalesfase=size(LFPtheta,2);
tquitar=0.5; 

wave1=nan(length(frecs1),Nptos,Nsenalesfase);
wave2=nan(length(frecs2),Nptos,Nsenalesamp);
for isen=1:size(LFPtheta,2)
    wave1(:,:,isen) = wavelet_mod(LFPtheta(:,isen),1/Fs,false,frecs1,6);
end
for isen=1:size(signals,2)
    wave2(:,:,isen) = wavelet_mod(signals(:,isen),1/Fs,false,frecs2,6);
end
fase1=angle(wave1(:,round(tquitar*Fs)+1:end-round(tquitar*Fs),:));
amp2=abs(wave2(:,round(tquitar*Fs)+1:end-round(tquitar*Fs),:));

Nbins=18; %tort
binedges=linspace(-pi,pi,Nbins+1); binedges(1)=-pi-0.1; binedges(Nbins+1)=pi+0.1;
tort=zeros(length(frecs1),length(frecs2),Nsenalesamp);
for if1=1:length(frecs1)  
    idbins=false(size(fase1,2),Nbins);
    for ibin=1:Nbins
        idbins(:,ibin)=and(fase1(if1,:)>binedges(ibin),fase1(if1,:)<(binedges(ibin+1)));
    end
    for if2=1:length(frecs2) 
        for isen=1:Nsenalesamp
            pamp=zeros(1,Nbins);
            for ibin=1:Nbins
                pamp(ibin)=mean(amp2(if2,idbins(:,ibin),isen));
            end
            pamp=pamp/sum(pamp);
            tort(if1,if2,isen)=(log(Nbins)+sum(pamp.*log(pamp)))/log(Nbins);
         end
    end
end

%%
fc1v=abs(frecs1-5);
fc1= find (fc1v==min(fc1v)); 
fc2v=abs(frecs1-15);  
fc2= find (fc2v==min(fc2v));
%fc1 = 4; fc2 = length(frecs1); % changed
fc3v=abs(frecs2-band(1));
fc3= find (fc3v==min(fc3v)); 
fc4v=abs(frecs2-band(2));
fc4= find (fc4v==min(fc4v)); 
tortOK=tort(fc1:fc2,fc3:fc4,:);
frecs2OK=frecs2(fc3:fc4);
% MI 
MI=zeros(Nsenalesamp,1);
    for i=1:Nsenalesamp
        MI(i,1)= mean2 (tortOK(:,:,i));
    end
% preffered freq gamma 
freq=zeros(Nsenalesamp,1);
aux3=zeros(length(tortOK(:,1,1,1)),Nsenalesamp);
for i=1:Nsenalesamp
    for j=1:length(tortOK(:,1,1,1))
        aux(j,i)=max(tortOK(j,:,i));
        aux2(j,i)=find(tortOK(j,:,i)==aux(j,i));
        aux3(j,i)=frecs2OK(aux2(j,i));
    end
    freq(i,1)=mean(aux3(:,i));
end

%%
% figure
% for isen=1:Nsenalesamp
% dibujar=abs(tort(:,:,isen));
% subplot(1,Nsenalesamp,isen); hold on;
% %surf(frecs1,frecs2,zeros(size(dibujar')), dibujar', 'EdgeColor', 'none', 'FaceColor', 'flat');
% contourf(frecs1,frecs2,dibujar',20);shading('flat');
% colorbar ('SouthOutside'); 
% %title(signals)
% xlim([frecs1(1) frecs1(end)])
% ylim([frecs2(1) frecs2(end)])
% if isen>1
%     set(gca,'YTick',[]);
% end
% end

