function X = gaussFilt(X,N1,varargin)

% Xsm = gaussFiltAng(X,N1)
% 
% 1 or 2-d gaussian filtering with border correction
% Input:
%     X: vector (or matrix) of angular values to smooth
%     N1: sd of gausswin for the first dimension. If N1==0, no smoothing
% Output:
%     Xsm: smoothed matrix
%     
% Adrien Peyrache 2011

if length(size(X))>2
    error('Too many dimensions!')
end

N2=0;
normBool = 1;
if ~isempty(varargin)
    N2 = varargin{1};
    if length(varargin)==2
        normBool = varargin{2};
    end
else
    N2=N1;
end

transpX = 0;
if size(X,1)==1
    X = X';
    transpX = 1;
end
   
l1 = size(X,1);

X = [X;X;X];
clear Xf

%s.d. is N/(2*alpha) and we want s.d. = original N.
%so here, with alpha=0.05, we multiply N by 10

    N1 = round(10*N1);
    gw = gausswin(N1,5); %alpha = 0.05
    if normBool
        gw = gw/sum(gw(:));
    end
    X = convn(X,gw,'same');

X = X(l1+1:2*l1,:);

if transpX
    X = X';
end
