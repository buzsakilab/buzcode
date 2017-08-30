function [peaktimes,winints,peakmags] = FindPeakInWin( data,timestamps,peakthresh,winthresh,minwindur,minseparation )
%[peaktimes,winints,peakmags] = FindPeakInWin(data,timestamps,peakthresh,winthresh,minwindur,minseparation )
%
%
%DLevenstein 2017
%% DEV
% data = -normGAMMA;
% timestamps = gammaLFP.timestamps;
% peakthresh = 1;
% winthresh = 0.5;
% minwindur = 0.03;
% minseparation = 0.01;

%%
[peakmags,peaktimes] = findpeaks(data,timestamps,'MinPeakHeight',peakthresh);

%%
overthresh = data>=winthresh;
crossup = find(diff(overthresh)==1);
crossdown = find(diff(overthresh)==-1);

%If there's only one crossing
if isempty(crossup) || isempty(crossdown)
    display('No windows over threshold')
    peaktimes = [];peakmags=[];winints=[];return
end

%Make sure first cross is up and last cross is down
if crossdown(1)<crossup(1); crossdown(1) = []; end
if crossup(end)>crossdown(end); crossup(end) = []; end

winints = timestamps([crossup crossdown]);

%Join and get rid of ints based on durations
winints = MergeSeparatedInts(winints,minseparation);
windur = diff(winints,[],2);
winints(windur<minwindur,:)=[];

[peaktimes,keepints,containingints,inint] = RestrictInts(peaktimes,winints);
winints = winints(containingints,:);
peakmags = peakmags(keepints);

%Get the biggest peak in each window
containingints = find(containingints);
for ii = 1:length(containingints)
    intpeaks = find(inint==containingints(ii));
    [biggestpeak,bigpeakidx] = max(peakmags(intpeaks));
    peakmags(ii) = biggestpeak;
    peaktimes(ii) = peaktimes(intpeaks(bigpeakidx));
end
peakmags((ii+1):end) = [];
peaktimes((ii+1):end) = [];

end

