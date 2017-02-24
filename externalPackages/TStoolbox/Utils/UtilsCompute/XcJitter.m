function [H0,Hm,HeI,HeS,Hstd,B] = XcJitter(t1,t2,binS,nbBins,varargin)

%  [H0,Hm,HeI,HeS,Hstd,B] = XcJitter(t1,t2,binS,nbBins,'optionname',optionvalue)
%  Computes jiterred cross-corr and returns an average and a intervalof
%  confidence
%  
%  INPUT:
%    t1,t2 : spike times in 10^-4 s
%    binS: bins ~size (in ms)
%    nbBins: number of bins
%
%   Options:
%    'confInt': width of the confidence interval (between 0 and 1, 0.95 default)
%    'nbIter': number of iterations (default 100)
%    'jitter': half width of jittering in ms (default 10 ms). Spk train
%    will be jittered in +/- this value homogeneously
%     
%  OUTPUT:
%    H0: normal CrossCorr
%    Hm: average jittered CrossCorr
%    HeI: lower bound of the interval of confidence
%    HeS: upper bound of the interval of confidence
%    Hstd = jittered cross-corr s.d. - (H0-Hm)./Hstd is thus the z-scored
%    Cross-Corr)
%    B: bin center points
%    
%  Adrien Peyrache 2010,2013 (big update)

% Parameters:
nbIter = 100;
confInt = 0.95;
jitter = 10; %in ms

for ii = 1:2:length(varargin),
  if ~isa(varargin{ii},'char'),
    error(['Parameter ' num2str(ii+4) ' is not a property (type ''help XcJitter'' for details).']);
  end
  switch(lower(varargin{ii})),
    case 'confint',
      confInt = varargin{ii+1};
      if ~isa(confInt,'numeric') | length(confInt) ~= 1 | confInt < 0,
        error('Incorrect value for property ''confInt'' (type ''help XcJitter'' for details).');
      end
    case 'nbiter',
      nbIter = varargin{ii+1};
      if ~isa(nbIter,'numeric') | length(nbIter) ~= 1 | nbIter < 0,
        error('Incorrect value for property ''nbIter'' (type ''help XcJitter'' for details).');
    end
    case 'jitter',
      jitter = varargin{ii+1};
      if ~isa(jitter,'numeric') | length(jitter) ~= 1 | jitter < 0,
        error('Incorrect value for property ''jitter'' (type ''help XcJitter'' for details).');
    end  
  end
end

if ~isempty(t1) && ~isempty(t2)
  [H0,B] = CrossCorr(t1,t2,binS,nbBins);
else
  error('Pb with vector size')
end

H1 = [];

for ii=1:nbIter
  %Generates a jittered spike train
  t2j = t2 + 2*10*jitter*(rand(length(t2),1)-0.5);
  t2j = sort(t2j);
  if ~isempty(t1) && ~isempty(t2j)
    H = CrossCorr(t1,t2j,binS,nbBins);
    H1 = [H1;H'];
  else
      display('Pb with vector size')
  end
end
Hm = mean(H1)';
HeI = zeros(length(B),1);
HeS = zeros(length(B),1);
Hstd = zeros(length(B),1);

for b=1:length(B)

  t = H1(:,b);
  t = sort(t,'ascend');
  HeI(b) = t(round((1-confInt)/2*nbIter));
  HeS(b) = t(round((confInt + (1-confInt)/2)*nbIter));
  Hstd(b) = std(t);
end


