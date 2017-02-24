function [SM, EM] = NonRippleTimes(S, E)
% [SM, EM] = RipplesMidPoint(S, E)
% 
% given the ripples times, gives time intervals placed randomly
% between ripple events, for sake of controls. The duration of each "mid"
% event is equal to the duration of the immediately previous ripple event.
% INPUTS:
% S: array with the ripples events start times
% E: array with the ripples events end times
%
% OUTPUTS:
% SM: array with the non-ripple events start times
% SM: array with the non-ripple events end times

% batta 2000
% status beta


dr = (E(1:end-1) - S(1:end-1))/2;
mp = E(1:end-1) + dr + (S(2:end) - E(1:end-1) - dr) .* rand(1,length(dr)); 
SM = mp - dr;
EM = mp + dr;