function [] = chirpDetector(audioFile)

[b a]=butter(3,[25000./150000 75000./150000],'bandpass'); % create filter for chirp range
c=1;

for i = 1:150000:length(y)
   spectrogram(double(y(i:i+150000)),1000,[],[],Fs,'yaxis');
   pause(.1)
   
   meanPower(c) = mean(filtfilt(b,a,double(y(i:i+150000))));
   stdPower(c) = std(filtfilt(b,a,double(y(i:i+150000))));
   c=c+1;
end

end