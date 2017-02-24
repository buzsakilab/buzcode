function Xsm = gausssmooth(X,N1,varargin)

% Xsm = gausssmooth(X,N1,N2)
% 
% Input:
%     X: vector (or matrix) to smooth
%     N1: sd of gausswin for the first dimension. If N1==0, no smoothing
%     N2: in case of X is a matrix, sd of smoothing on the 2nd dimension. If omitted, N2=N1
%     
% Output:
%     Xsm: smoothed matrix
%     
% Adrien Peyrache 2011

warning('Now gausssmooth works in s.d., not number of samples')
if N1~=0
    N1 = 5*N1; %assuming alpha==2
else
    N1 = 1;
end
gw = gausswin(N1); 

if length(size(X))>2
    error('Too many dimensions!')
end

if min(size(X))>1
    if ~isempty(varargin)
        N2 = varargin{1};
        if N1~=0
            N2 = 5*N2; %assuming alpha==2
        else
            N2 = 1;
        end
    else
        N2=N1;
    end
    gw = gw*gausswin(N2)';
else
    X = X(:);
end

gw = gw/sum(gw(:));
Xsm = convn(X,gw,'same');
 