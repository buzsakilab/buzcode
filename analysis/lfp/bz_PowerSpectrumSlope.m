function [specslope,spec] = bz_PowerSpectrumSlope(lfp,dt,winsize)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Gao, Peterson, Voytek 2016

%%
% winsize = 4; 
% dt = 0.5;
%%
%Calcluate spectrogram
noverlap = winsize-dt;
spec.freqs = logspace(0.5,2,200);
winsize = winsize .*lfp.samplingRate;
noverlap = noverlap.*lfp.samplingRate;
[spec.data,~,spec.timestamps] = spectrogram(lfp.data,winsize,noverlap,spec.freqs,lfp.samplingRate);
spec.amp = log10(abs(spec.data));

%% Fit the slope of the power spectrogram
rsq = zeros(size(spec.timestamps));
s = zeros(length(spec.timestamps),2);
yresid = zeros(length(spec.timestamps),length(spec.freqs));
for tt = 1:length(spec.timestamps)
    %Fit the line
    x = log10(spec.freqs);  y=spec.amp(:,tt)';
    s(tt,:) = polyfit(x,y,1);
    %Calculate the residuals
    yfit =  s(tt,1) * x + s(tt,2);
    yresid(tt,:) = y - yfit;
    %Calculate the rsquared value
    SSresid = sum(yresid(tt,:).^2);
    SStotal = (length(y)-1) * var(y);
    rsq(tt) = 1 - SSresid/SStotal;
end

%% Output Structure
specslope.data = s(:,1);
specslope.intercept = s(:,2);
specslope.timestamps = spec.timestamps;
specslope.samplingRate = 1./dt;
specslope.winsize = winsize;

specslope.rsq = rsq;
specslope.resid = yresid;

end

