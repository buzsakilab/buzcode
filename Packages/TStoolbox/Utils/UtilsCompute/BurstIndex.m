function [burstBool burstIsiPk] = BurstIndex(S,varargin)

% [burstBool burstIsiPk] = BurstIndex(S,ep)
% 
% Computes burst index
% Input: 
%     S: a ts object of spike timing.
%     ep (optionnal): intervalSet object of the epoch of interest
%     
% Output:
%     burstIdx: burstiness index of the cell (1 for a Poisson process)

%Adrien Peyrache

burstLimUp = 8; %max isis in burst, in ms
burstLimDw = 1.5; %lower bound of isis burst mode, in ms

ratioBurst = 0.05; %minimum ratio of observed over expected isi histogram peak.
peakSeparation = 0.5; %in number of log-scale decades
db = 0.05; %bin size of ISI histogram
burstBool = 0;
burstIsiPk = NaN;


b = [-1:db:4];

if ~isempty(varargin)
    ep = varargin{1};
else
    ep = timeSpan(S);
end

s = Range(Restrict(S,ep),'ms');
dt = diff(s);
dt = dt(dt<10^b(end));
h = hist(log10(dt),b(1:end-1));
hs = gaussFilt(h(:),2,0)/(length(s)-1);
r = rate(S,ep);

r1 = r/1000;

expDt = r1*exp(-r1.*10.^b(1:end-1)).*diff(10.^b);
expDt = expDt(:);

ix  = LocalMinima(-hs,peakSeparation/db,0);
ix(b(ix)<log10(burstLimDw)) = [];

if ~isempty(ix)
%     disp(10^b(ix(1)) < burstLimUp)
%     disp(hs(ix(1))./expDt(ix(1))>ratioBurst)
    n = sum(dt>=burstLimDw & dt<=burstLimUp)/length(dt);
    burstBool = 10^b(ix(1)) < burstLimUp & n>ratioBurst;
    burstIsiPk = 10^b(ix(1));
    
end