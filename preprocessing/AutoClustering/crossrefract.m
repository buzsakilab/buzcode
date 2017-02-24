function [fxy fcorr] = crossrefract(rgx,rgy,tR,tC);

% USAGE:
%   [fxy fcorr] = CROSSREFRACT(rgx,rgy,tR,tC)
%   
% INPUTS:
%   
% CROSSREFRACT computes the fraction fxy of rogue spike in the merged
% spike trains. fcorr is an attempt to normalize it relatively to the
% number of spikes in each individual spike train (not used so far)
% 
% Adrien Peyrache, 2012

lx = length(rgx);
fx = FractionRogueSpk(rgx,tR,tC);
ly = length(rgy);
fy = FractionRogueSpk(rgy,tR,tC);
fxy = FractionRogueSpk(sort([rgx;rgy]),tR,tC);
fcorr = fxy-(lx*fx+ly*fy)/(lx+ly);
