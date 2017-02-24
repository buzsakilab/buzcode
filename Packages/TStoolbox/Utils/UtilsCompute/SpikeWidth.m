function [p2v, hp, hv, pkratio] = SpikeWidth(meanWaveF,Fs)

% USAGE:
%     [peak2val,halfPkWidth] = SpikeWidth(meanWaveF,Fs)
% 
% INPUT:
%     meanWaveF: a N by M matrix of a single cell averaved waveforms from N
%     electrodes , M samples
%     Fs: sampling rate (in Hz)
%     

nbE = size(meanWaveF,1);
if nbE>1
    spkPow = sum(meanWaveF.^2,2);
    [dummy,mxIx] = max(spkPow);
    w = meanWaveF(mxIx,:);
else
    w = meanWaveF;
end

hp=0;
hv=0;
p2v=0;
pkratio=0;

if ~any(isnan(w))
        wu = w;
    	wu = resample(wu,10,1);
        nSamples = length(wu);

    	t = [0:1/(Fs*10):(nSamples-1)/(Fs*10)]*1000;
        [minVal,minPos] = min(wu);
    	maxPos1 = LocalMinima(-wu(1:minPos),10,0);
       	maxPos2 = LocalMinima(-wu(minPos+1:end),10,0);
        
        if ~isempty(maxPos1) && ~isempty(maxPos2)
        
        maxPos1 = maxPos1(end);
        maxPos2 = maxPos2(1);
        maxPos2 = maxPos2+minPos;
        
        maxVal1 = wu(maxPos1);
        maxVal2 = wu(maxPos2);
        
        pkratio = (maxPos2-maxPos1)/(maxPos2+maxPos1);
        
    	wun = wu/maxVal2;
        ix = match(0.5,wun(minPos+1:maxPos2),0.05,1);
        wun = wu/minVal;
        ix1 = match(0.5,wun(maxPos1+1:minPos),0.05,1);
        ix2 = match(0.5,wun(minPos+1:maxPos2),0.05,1);
        
    	if ~isempty(ix)
          hp = 2*(t(maxPos2)-t(maxPos2-ix));
          hv = t(minPos+ix2)-t(maxPos1+ix1);
    	  p2v = t(maxPos2)-t(minPos);
          pkratio = (maxPos2-maxPos1)/(maxPos2+maxPos1);
        end
        end
end
  