function [xy, c_lo, c_hi] = nancorrcoef(varargin)
%[xy, clo, chi] = nancorrcoef(x,y, alpha) Correlation coefficients.
%   NANCORRCOEF(X) is a matrix of correlation coefficients formed
%   from array X whose each row is an observation, and each
%   column is a variable.
%   NANCORRCOEF(X,Y), where X and Y are column vectors is the same as
%   NANCORRCOEF([X Y]).
%   
%   If C is the covariance matrix, C = COV(X), then CORRCOEF(X) is
%   the matrix whose (i,j)'th element is
%
%          C(i,j)/SQRT(C(i,i)*C(j,j)).
%
% correlation matrix is computed from only those datapoints that had all
% finite data
% INPUTS: 
% X, or X, Y, the data
% alpha: the confidence level for confidence interval
% ci_method: the method used to compute the confidence interval, can be
% "fisher" or "bootstrap"  
% OUTPUTS:
% xy: the correlation coefficient
% Clo, Chi (optional): the confidence interval for r, computed by means of
% Fisher transform (default; see Sheskin p. 779) or by bootstrap resampling

% Reference
% Sheskin D.J. 2000 Handbook of parametric and non-parametrix statistical
% procedures, 2nd Ed. Chapman & Hall
  
  
  
ci_method = 'bootstrap';
alpha = 0.05;

narg = nargin;


if ischar(varargin{end})
  ci_method = varargin{end};
  narg = narg-1;
end






if narg == 0
  error('not enough input arguments');
end


if narg >= 1
  x = varargin{1};
end

if narg >= 2
  y = varargin{2};
end

if narg >= 3
  alp = varargin{3};
end



if narg == 1
  
  D = x;
end

if narg == 2 
  
  if  (size(x, 1) == size(y, 1))
    D = [x, y];
  else
    D = x;
    alpha = y;
    end
end


if narg == 3
  D = [x, y];
  alpha = alp;
end





ix = find(isfinite(sum(D, 2)));

D = D(ix, :);
c = cov(D);

d = diag(c);
xy = c./sqrt(d*d');



if nargout == 3
  
  switch ci_method
   case 'fisher'
    n = size(D, 1);
    z = atanh(xy);
    zc = erfinv(1 - alpha) * sqrt(2/(n-3));
    c_hi = tanh(z +zc);
    c_lo = tanh(z -zc);
    case 'bootstrap'
    n = size(D,1);
    dim = size(D, 2);
    nboot = 100/alpha;
    boot = zeros(size(D,2), size(D,2), nboot);
    for i = 1:nboot
      bt_idx = ceil(n * rand(n,1));
      boot(:,:,i) = corrcoef(D(bt_idx,:));
    end
    
    c_hi = diag(ones(dim,1));
    c_lo = diag(ones(dim,1));
    for i = 1:dim-1
      for j = (i+1):dim
	bx = sort(squeeze(boot(i,j,:)));
	c_hi(i,j) = bx(floor(nboot*(1-(alpha/2))));
	c_hi(j,i) = c_hi(i,j);
	c_lo(i,j) = bx(floor(nboot*(alpha/2)));
	c_lo(j,i) = c_lo(i,j);
      end
    end
    
   
   otherwise
    error('unknown Ci method');
  end
end

    
    
    
    
  