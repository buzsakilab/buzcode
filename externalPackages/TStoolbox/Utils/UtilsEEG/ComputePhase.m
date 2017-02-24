function thetaPhase = computePhase(T, tPeaks)
% compute theta phase of spikes 
% 

% copyright (c) 1999-2006 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html



lcell = Data(T);
  
  for i = 1:length(lcell)
    tp = binsearch_floor(tPeaks, lcell(i));
    thetaPhase(i) = (lcell(i) -  tPeaks(tp)) / (tPeaks(tp+1) - ...
	tPeaks(tp));
    if thetaPhase(i) >= 1 | thetaPhase(i) < 0
      error('Phase out of range');
    end
  end