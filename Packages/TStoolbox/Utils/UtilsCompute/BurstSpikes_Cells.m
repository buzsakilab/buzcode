function [spkBurstIx,burstIdx,burstLg] = BurstSpikes(S,ep)

Burst = cell(length(S),1);
burstCell = zeros(length(S),1);

burstLim = 10; %in ms

S = Restrict(S,ep);
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

expectedBurstIdx = (1-exp(-(burstLim/1000)*rate(S,ep)));

% The final burstIdx is then:
burstIdx = max(sum(spkBurstIx>0)/length(rg) - expectedBurstIdx,0);
