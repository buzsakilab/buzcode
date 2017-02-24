function [SM, EM] = RipplesMidPoint(S, E)
% [SM, EM] = RipplesMidPoint(S, E)
% 
% given the ripples times, gives time intervals placed at the mid-points
% between ripple events, for sake of controls. The duration of each "mid"
% event is equal to the duration of the immediately previous ripple event.

mp = (S(2:end) + E(1:end-1))/ 2;
dr = (E(1:end-1) - S(1:end-1))/2;

SM = mp - dr;
EM = mp + dr;