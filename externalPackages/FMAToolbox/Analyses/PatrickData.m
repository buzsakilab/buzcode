function PatrickData(CNO4lfp,SAL4lfp)

figure;
for i=1:4
  subplot(4,1,i);
  MTSpectrogram([CNO4lfp(:,1) CNO4lfp(:,i+1)],'show','on','range',[0 80]);
  ylabel('Freq');
  clim([0 10]);
end
xlabel('Time');
suptitle('CNO - amy1/cx1/amy2/cx2');

figure;
for i=1:4
  subplot(4,1,i);
  MTSpectrogram([SAL4lfp(:,1) SAL4lfp(:,i+1)],'show','on','range',[0 80]);
  ylabel('Freq');
  clim([0 10]);
end
xlabel('Time');
suptitle('Saline - amy1/cx1/amy2/cx2');

for i=[2 4]
  MTCoherogram([CNO4lfp(:,1) CNO4lfp(:,i)],[CNO4lfp(:,1) CNO4lfp(:,i+1)],'show','on','range',[0 40]);
  xlabel(['Coherence - CNO ' int2str(i)]);
end

for i=[2 4]
  MTCoherogram([SAL4lfp(:,1) SAL4lfp(:,i)],[SAL4lfp(:,1) SAL4lfp(:,i+1)],'show','on','range',[0 40]);
  xlabel(['Coherence - SAL ' int2str(i)]);
end

nfft= 2^nextpow2(2*4000); %Number of FFT points to be used
fspec = linspace(0, 4000/2, nfft/2 + 1); %Frequency Vector
for i=[2 4]
  cohspectrum = mscohere(CNO4lfp(:,i),CNO4lfp(:,i+1),hamming(nfft), 0.5*nfft, nfft, 4000);
  figure;
  plot(fspec(fspec<100),cohspectrum(fspec<100));
  xlabel('Coherence Spectrum - CNO');
end
for i=[2 4]
  cohspectrum = mscohere(SAL4lfp(:,i),SAL4lfp(:,i+1),hamming(nfft), 0.5*nfft, nfft, 4000);
  figure;
  plot(fspec(fspec<100),cohspectrum(fspec<100));
  xlabel('Coherence Spectrum - SAL');
end



