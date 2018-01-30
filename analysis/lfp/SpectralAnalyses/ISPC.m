function [ IPSC,phaselag ] = ISPC(sig1phase,sig2phase)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
%   sig2phase can be the phase of signal 2 or 'diff', indicating that
%           sig1phase is the phase difference between two signals
%
%%

if strcmp(sig2phase,'diff')
    phasediff = sig1phase;
else
    phasediff = sig1phase-sig2phase;
end

phasesynch = (mean(exp(1i.*phasediff),1));
IPSC = abs(phasesynch);
phaselag = angle(phasesynch);



end

