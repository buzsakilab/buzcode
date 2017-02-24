function [meanwav, stdwav, meansw, stdsw, meanp2v, stdp2v, peakheight] = WaveFormsStats(TTfilename, T)

chunk_length = 5000;
offset = 1;
meanwav = zeros(4,32);
stdwav = zeros(4,32);
meansw = zeros(1,4);
stdsw = zeros(1,4);
meanp2v = zeros(1,4);
stdp2v = zeros(1,4);
peakheight = zeros(0,4);
nSpikes = 0;


while 1
  TT = LoadTT(TTfilename, offset, chunk_length);
  nReadSpikes = length(Range(TT,'ts'));
  if nReadSpikes == 0
    break
  end;
  EndingTime = EndTime(TT);
  StartingTime = StartTime(TT);
  T1 = Restrict(T, StartingTime, EndingTime);
  TT = Restrict(TT, Data(T1));
  
  
  
  
  t = Range(TT, 'ts');
  wv = Data(TT);
  nReadSpikes = length(t);

  for i = 1:4
    meanwav(i,:) = meanwav(i,:) + squeeze(sum(wv(:,i,:), 1))';
    stdwav(i,:) = stdwav(i,:) + squeeze(sum(wv(:,i,:).*wv(:,i,:), 1))';
  end

  % spike width and peaktovalley
  [peak, ipeak] = max(wv, [], 3);
  [vlly, ivlly] = min(wv, [], 3);
  peakheight = cat(1, peakheight, peak);
  sw = ivlly - ipeak;
  meansw = meansw +sum(sw, 1);
  stdsw = stdsw + sum(sw .*sw, 1);
  S = warning;
  warning off;
  p2v = abs(peak)./abs(vlly);
  warning(S);
  p2v(find(vlly==0)) = 0;
  meanp2v = meanp2v+sum(p2v, 1);
  stdp2v = stdp2v +sum(p2v,1);
  
  nSpikes = nSpikes + nReadSpikes;
  offset = offset + chunk_length;
end

nSpikes

meanwav = meanwav / nSpikes;
meansw = meansw / nSpikes;
meanp2v = meanp2v / nSpikes;
stdwav = stdwav / (nSpikes-1);
stdsw = stdsw / (nSpikes-1);
stdp2v = stdp2v / (nSpikes-1);
stdwav = stdwav - (nSpikes/(nSpikes-1)) * meanwav .* meanwav;
stdsw = stdsw - (nSpikes/(nSpikes-1)) * meansw .* meansw;
stdp2v = stdp2v - (nSpikes/(nSpikes-1)) * meanp2v .* meanp2v;
stdwav = sqrt(stdwav);
stdsw = sqrt(stdsw);
stdp2v = sqrt(stdp2v);

