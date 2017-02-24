function [meanwav, stdwav] = MeanStdWaveforms(TTfilename, T)

chunk_length = 5000;
offset = 1;
meanwav = zeros(4,32);
stdwav = zeros(4,32);
nSpikes = 0;

while 1
  TT = LoadTT(TTfilename, offset, chunk_length);
  EndingTime = EndTime(TT);
  StartingTime = StartTime(TT);
  T1 = Restrict(T, StartingTime, EndingTime);
  TT = Restrict(TT, Data(T1));
  
  
  
  
  t = Range(TT, 'ts');
  wv = Data(TT);
  nReadSpikes = length(t);

  
  if nReadSpikes == 0
    break
  end;
  for i = 1:4
    meanwav(i,:) = meanwav(i,:) + squeeze(sum(wv(:,i,:), 1))';

  end

  
  
  nSpikes = nSpikes + nReadSpikes;
  offset = offset + chunk_length;
end

nSpikes
meanwav = meanwav / nSpikes;

offset = 1;
nSpikes = 0;
while 1
  TT = LoadTT(TTfilename, offset, chunk_length);
  EndingTime = EndTime(TT);
  StartingTime = StartTime(TT);
  T1 = Restrict(T, StartingTime, EndingTime);
  TT = Restrict(TT, Data(T1));
  t = Range(TT, 'ts');
  wv = Data(TT);
  nReadSpikes = length(t);

  
  if nReadSpikes == 0
    break
  end
  
  m = repmat(meanwav, [1 1 nReadSpikes]);
  m = permute(m, [3 2 1]);
  
  wv = wv-1;
  
  
  for i = 1:4
    stdwav(i,:) = stdwav(i,:) + squeeze(sum(wv(:,i,:).*wv(:,i,:), 1))';

  end

  
  
  nSpikes = nSpikes + nReadSpikes;
  offset = offset + chunk_length;
end


stdwav = stdwav / (nSpikes-1);
stdwav = sqrt(stdwav);

