%% determine STFT parameters

hop = 128;
nfft = 1024;
win = hann(nfft);

%% load signals

load('handel.mat','y','Fs')
xt = y;
load('laughter.mat','y','Fs')
xi = y;

% crop xt
xt = xt(1:length(xi));

%% calculate and plot STFTs

% calculate STFTs
[st,f,t] = iosr.dsp.stft(xt,win,hop,Fs);
si = iosr.dsp.stft(xi,win,hop,Fs);

% calculate spectrograms
stlog = 20.*log10(abs(st));
silog = 20.*log10(abs(si));

% calculate colormap limits
clim = [min([stlog(:); silog(:)]) max([stlog(:); silog(:)])];

% draw figures
figure
%
subplot(8,4,[1 2 5 6])
imagesc(t,f,stlog)
set(gca,'ydir','normal','xticklabel',[],'clim',clim)
title('Target spectrogram')
ylabel('Frequency [Hz]')
ch = colorbar;
ylabel(ch,'Magnitude [dB]')
%
subplot(8,4,[9 10 13 14])
imagesc(t,f,silog)
set(gca,'ydir','normal','clim',clim)
title('Interference spectrogram')
xlabel('Time [s]')
ylabel('Frequency [Hz]')
ch = colorbar;
ylabel(ch,'Magnitude [dB]')

%% calculate and plot ideal masks

% calculate ideal masks
[irm,ibm] = iosr.bss.idealMasks(st,si);

% draw figures
subplot(8,4,[3 4 7 8])
imagesc(t,f,irm)
set(gca,'ydir','normal','xticklabel',[],'yticklabel',[],'clim',[0 1])
title('Ideal ratio mask (IRM)')
ch = colorbar;
ylabel(ch,'Mask value')
%
subplot(8,4,[11 12 15 16])
imagesc(t,f,ibm)
set(gca,'ydir','normal','yticklabel',[],'clim',[0 1])
title('Ideal binary mask (IBM)')
xlabel('Time [s]')
ch = colorbar;
ylabel(ch,'Mask value')

%% resynthesise and plot output

% calculate and apply ideal masks in one step
[z_irm,z_ibm,t] = iosr.bss.applyIdealMasks(xt,xi,win,hop,Fs);

% draw figures
subplot(8,4,21:24)
plot(t,xt)
set(gca,'xticklabel',[],'xlim',[0 max(t)],'ylim',[-1 1])
title('Clean target signal')
ylabel('Amplitude')
%
subplot(8,4,25:28)
plot(t,z_irm)
set(gca,'xticklabel',[],'xlim',[0 max(t)],'ylim',[-1 1])
title('Target resynthesised from IRM')
ylabel('Amplitude')
%
subplot(8,4,29:32)
plot(t,z_ibm)
set(gca,'xlim',[0 max(t)],'ylim',[-1 1])
title('Target resynthesised from IBM')
xlabel('Time [s]')
ylabel('Amplitude')
