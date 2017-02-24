function [S, E] = findTheta(Filt_EEG, HighThreshold, LowThreshold, varargin)
% [S, E] = findTheta(Filt_EEG, HighThreshold, LowThreshold) find theta events
% 
% INPUTS: 
% Filt_EEG: EEG tsd filtered in the theta range
% HighThreshold: threshold for detection of transitions from LIA to THETA
% LowThreshold: Threshold for detection of transitions from THETA to LIA
%
% PARAMETERS:
% Qup: number of cycles to transit to THETA
% Qdown: number of cycles to transit to LIA
% CloseThreshold: Closeness threhshold,if two theta events are closer than
%   this,  lump them together
% MinThetaDuration: discard theta evants shorter than this
% OUTPUTS:
% S: ts object with theta events start times
% E: ts object with theta events end times

% batta 1999
% status: beta



Qup = 10;
Qdown = 5;
CloseThreshold = 20000 * 10; % if two theta events are closer than this, lump
                            % them together



MinThetaDuration  = 10000* 5; % discard theta evants shorter than this

%  Extract_varargin;

LIA = -1;
THETA = 1;

state = 0;
% do it in chunks;

EEGStart = StartTime(Filt_EEG);
EEGEnd = EndTime(Filt_EEG);

% chunk length in timestamps

LChunk = 3000000;

ChStart = EEGStart;
ChEnd = min(ChStart + LChunk, EEGEnd);
nChunks = (EEGEnd - EEGStart)/LChunk;
TTStart = [];
TTEnd = [];

ii = 0;
%  h = waitbar(0, 'Find Ripples...');

while ChEnd <= EEGEnd

%    waitbar(ii/nChunks, h);

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
  PeaksIdx = find(de1 .* de2 < 0);
  
  Peaks = eeg(PeaksIdx);
  Peaks = abs(Peaks);
  
  if state == 0
    if Peaks(1:3) > HighThreshold
      state = THETA;
    else
      state = LIA;
    end
  end
  
  % when  Qup peaks in a row are greater than HighThreshold we detect an
  % "upwards" transition, that is a transition from LIA to THETA
  
  % when Qdown peaks in a row are smaller than LowThreshold, we detect a
  % "downwards" transition, from THETA to LIA
  
  q = max(Qup, Qdown);
  P = zeros(q, length(Peaks)+q-1);
  
  P(1,1:length(Peaks)) = Peaks;
  for i = 2:q
    P(i,2:end) = P(i-1,1:end-1);
  end
  
  Pmax = max(P(1:Qdown, :));
  Pmin = min(P(1:Qup,:));
  
  Pmax = Pmax(Qdown-1:end);
  Pmax = Pmax(1:length(Peaks));
  Pmin = Pmin(Qup-1:end);
  Pmin = Pmin(1:length(Peaks));
  
  UpTransIdx = find(Pmin > HighThreshold);
  UpTransIdx = PeaksIdx(UpTransIdx);
  DownTransIdx = find(Pmax < LowThreshold);
  DownTransIdx = PeaksIdx(DownTransIdx);
  
  % then we define a theta period as the period beteween an "Upwards"
  % transition and the next "Downwards" transition
  
  Start_i = 1;
  End_i = 1;
  StartThetaIdx = [];
  EndThetaIdx = [];
  CP = 0;
  while 1
    if state == LIA
      U = UpTransIdx(find((UpTransIdx > CP)));
      if isempty(U)
	break;
      end
      
      StartThetaIdx(Start_i) = min(U);
      CP = StartThetaIdx(Start_i);
      Start_i = Start_i +1;
      state = THETA;
    else
      D = DownTransIdx(find(DownTransIdx > CP));
      if isempty(D)
	break;
      end
      
      EndThetaIdx(End_i) = min(D);
      CP = EndThetaIdx(End_i);
      End_i = End_i +1;
      state = LIA;
    end
  end  

  
  TTStart = [TTStart t(StartThetaIdx)];
  TTEnd = [TTEnd t(EndThetaIdx)];

  ii = ii+1;
  ChStart = ChStart + LChunk;
  ChEnd = ChEnd + LChunk;

end
%  close(h);
i = 2;


if state == THETA
  TTEnd = [TTEnd EEGEnd];
end




if TTEnd(1) <TTStart(1)
  TTStart = [EEGStart TTStart];
  TTEnd = [TTEnd EEGEnd];
end


% lumping together theta events closer than CloseThreshold

while i <= length(TTStart)
  if (TTStart(i) - TTEnd(i-1)) < CloseThreshold
    TTStart = [TTStart(1:i-1) TTStart(i+1:end)];
    TTEnd = [TTEnd(1:i-2) TTEnd(i:end)];
  else
    i = i +1 ;
  end
  
end

% discarding theta events shorter that MinThetaDuration

i = 1;
while i <= length(TTStart)
  if(TTEnd(i) - TTStart(i)) < MinThetaDuration
    TTStart = [TTStart(1:i-1) TTStart(i+1:end)];
    TTEnd = [TTEnd(1:i-1) TTEnd(i+1:end)];
  else
    i = i+ 1;
  end
  
end

  
S = ts(TTStart);
E = ts(TTEnd);