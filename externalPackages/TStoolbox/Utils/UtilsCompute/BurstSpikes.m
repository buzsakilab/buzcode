function [spkBurstIx,burstLg,burstIdx] = BurstSpikes(S,varargin)
% 
% [spkBurstIx,burstLg] = BurstSpikes(S)
% 
% Detects burst spikes and computes burst index
% Input: 
%     S: a ts object of spike timing.
% Output:
%     spkBurstIx: a vector same length as S, indicating spike index in each
%     burst (0 for isolated spikes)
%     burstLg: number of spikes per burst
%
% [spkBurstIx,burstLg] = BurstSpikes(S,burstISI)
%     (optional): determines the maximum ISI in a burst. Default is 10ms
%
% [spkBurstIx,burstLg,burstIdx] = BurstSpikes(S,burstISI,ep)
%     outputs the burstiness index of the cell (0 for a Poisson process)
%     The intervalSet 'ep' could be put as an additional argument to compute
%     the firing rate. Otherwise, it is computed on the interval bounded by 
%     the first and last spikes.

%Adrien Peyrache

burstLim = 10; %in ms
if ~isempty(varargin)
    burstLim = varargin{1};
    if ~isnumeric(burstLim)
        error('Burst isi must be numeric')
    end
end
rg = Range(S,'ms');
dr = diff(rg);

shortDt = dr < burstLim;

spkBurstIx = zeros(length(rg),1,'int8');
burstLg = [];

u=1;
while u<length(shortDt)
   if shortDt(u)
        spkBurstIx(u) = 1;
        spkBurstIx(u+1) = 2;
        v=1;
        while shortDt(u+v) && (u+v)<length(shortDt)
          spkBurstIx(u+v+1) = v+2;
          v=v+1;
        end
        burstLg = [burstLg;v];
        u = u+v;
   end
   u = u+1;
end

%Athe expected burst index for a Poisson process is:
if nargout==3
    if length(varargin)==2
        ep = varargin{2};
        expectedBurstIdx = (1-exp(-(burstLim/1000)*rate(S,ep)));
    else
        expectedBurstIdx = (1-exp(-(burstLim/1000)*rate(S)));
    end
    
    % The final burstIdx is then:
    burstIdx = max(sum(spkBurstIx>0)/length(rg) - expectedBurstIdx,0);

end