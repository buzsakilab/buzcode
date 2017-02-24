function [meanwav, stdwav, meansw, stdsw, meanp2v, stdp2v, peakheight, GoodSpikes] = WaveFormsStats(TTfilename, T);

chunk_length = 20000;
offset = 1;
meanwav = zeros(4,32);
stdwav = zeros(4,32);
meansw = zeros(1,4);
stdsw = zeros(1,4);
meanp2v = zeros(1,4);
stdp2v = zeros(1,4);
peakheight = zeros(0,4);
nSpikes = 0;
StartingTime = 0;
nBadSpikes = 0;
GoodSpikes = [];

t = Data(T);
dt = diff(t);
dt = find(dt >  0) + 1;
t = t(dt);
T = ts(t);
while 1
  [TT, TTbad] = LoadTTbatta(TTfilename, offset, chunk_length);
  nReadSpikes = length(Range(TT,'ts'));
  if nReadSpikes == 0
    break
  end;
  EndingTime = EndTime(TT)-20;
  if EndingTime < StartingTime
    error('Problme with TT file segmentation');
  end
  
  StartingTime = StartTime(TT);
  T1 = Restrict(T, StartingTime, EndingTime);
  t1 = Data(T1);
  tbad = Range(TTbad, 'ts');
  z = zeros(size(t1));
  for l = 1:length(tbad)
    z(find(t1==tbad(l))) = 1;
  end
  nBadSpikes = nBadSpikes + sum(z);
  t1 = t1(find(z == 0));
    GoodSpikes = [GoodSpikes;t1];
  TT = Restrict_adr(TT, t1);
  
  
  
  
  t = Range(TT, 'ts');
  if length(find(t ~= t1)) ~= 0
    length(find(t ~= t1))
    error('Problem with Restrict');
  end
  
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
nBadSpikes
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


