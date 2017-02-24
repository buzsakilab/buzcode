function [S, E, M] = findRipples(Filt_EEG, DetectThreshold, LimitThreshold, varargin)
% [S, E, M] = findRipples(Filt_EEG, DetectThreshold, LimitThreshold, varargin)
% 
% INPUTS: 
% Filt_EEG: EEG tsd filtered in the ripples (e.g. 100-300 Hz) range
% DetectThreshold: Threshold for detection of ripples 
% LimitThreshold: Threshold for finding the ripple boundaries
%
% PARAMETERS:
% Q1: number of cycles to check to find boundaries
% CloseThreshold: Closeness threhshold,if two ripple events are closer than
%   this,  lump them together
% MinRippleDuration: discard theta evants shorter than this
% OUTPUTS:
% S: ts object with ripple events start times
% E: ts object with ripple events end times
% M: ts object with ripple events peak times

% batta 1999
% status: beta


% parameters
Q1 = 3;
CloseThreshold = 50 * 10;
MinRippleDuration = 30 * 10;

%  Extract_varargin;
% do it in chunks;
MinRippleDuration;

EEGStart = StartTime(Filt_EEG);
EEGEnd = EndTime(Filt_EEG);

% chunk length in timestamps

LChunk = 3000000;

ChStart = EEGStart;
ChEnd = min(ChStart + LChunk, EEGEnd);
nChunks = (EEGEnd - EEGStart)/LChunk;
TRStart = [];
TREnd = [];
TRMax = [];
TRMaxValue = [];

i = 0;
h = waitbar(0, 'Find Ripples...');

while ChStart < EEGEnd
  waitbar(i/nChunks, h);

  Ch = Restrict(Filt_EEG, ChStart, ChEnd);
  eeg = Data(Ch);
  t = Range(Ch, 'ts');
  sz = size(eeg);
  if sz(1) ~= 1
    eeg = eeg';    % we want to work with row vectors
    t = t';
  end
  
  

  
  de = diff(eeg);
  de1 = [de 0];
  de2 = [0 de];
  
  
  %finding peaks
  upPeaksIdx = find(de1 < 0 & de2 > 0);
  downPeaksIdx = find(de1 > 0 & de2 < 0);
  
  PeaksIdx = [upPeaksIdx downPeaksIdx];
  PeaksIdx = sort(PeaksIdx);
  
  Peaks = eeg(PeaksIdx);
  Peaks = abs(Peaks);
  
  %when a peaks is above threshold, we detect a ripple
  
  RippleDetectIdx = find(Peaks > DetectThreshold);
  DetectDiff = [0 diff(RippleDetectIdx)];
  RippleDetectIdx = RippleDetectIdx(find(DetectDiff > 2));
  RippleDetectIdx = RippleDetectIdx(find(RippleDetectIdx < length(Peaks)-Q1+1));
  RippleStart = zeros(1, length(RippleDetectIdx));
  RippleEnd = zeros(1, length(RippleDetectIdx));
  RippleMax = zeros(1, length(RippleDetectIdx));
  RippleMaxValue = zeros(1, length(RippleDetectIdx));
  for ii = 1:length(RippleDetectIdx) 
    CP = RippleDetectIdx(ii); % Current Peak
    % detect start of the ripple
    for j = CP-1:-1:Q1
      if all(Peaks(j-Q1+1:j) < LimitThreshold)
	
	break;
      end
    end
    RippleStart(ii) = j;
    %detect end of ripple
    for j = CP+1:length(Peaks)-Q1+1
      if all(Peaks(j:j+Q1-1)< LimitThreshold)
	
	break;
      end
    end
    RippleEnd(ii) = j;
    [RippleMaxValue(ii), RippleMax(ii)] = max(Peaks([RippleStart(ii):RippleEnd(ii)]));
    RippleMax(ii) = RippleStart(ii) + RippleMax(ii) - 1;
  end
  
  TRStart = [TRStart t(PeaksIdx(RippleStart))];
  TREnd = [TREnd t(PeaksIdx(RippleEnd))];
  TRMax = [TRMax t(PeaksIdx(RippleMax))];
  TRMaxValue = [TRMaxValue RippleMaxValue];
  i = i+1;
  ChStart = ChStart + LChunk;

  ChEnd = min(ChStart + LChunk, EEGEnd);
end
close(h);
i = 2;


while i <= length(TRStart)
  
  if (TRStart(i) - TREnd(i-1)) < CloseThreshold
    TRStart = [TRStart(1:i-1) TRStart(i+1:end)];
    TREnd = [TREnd(1:i-2) TREnd(i:end)];
    [mx, ix] = max([TRMaxValue(i-1) TRMaxValue(i)]);
    TRMax = [TRMax(1:i-2) TRMax(i - 2 + ix) TRMax(i+1:end)]; 
  else
    i = i +1 ;
    
  end
  
end

% control
length(TRStart)
length(TREnd)

i = 1;
while i <= length(TRStart)
  if(TREnd(i) - TRStart(i)) < MinRippleDuration
    
    TREnd(i) - TRStart(i);
    
    TRStart = [TRStart(1:i-1) TRStart(i+1:end)];
    TREnd = [TREnd(1:i-1) TREnd(i+1:end)];
    TRMax = [TRMax(1:i-1), TRMax(i+1:end)];
  else
    i = i+ 1;
    
  end
  
end

  
S = ts(TRStart);
E = ts(TREnd);
M = ts(TRMax);