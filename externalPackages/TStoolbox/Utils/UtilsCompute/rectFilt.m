function X = gaussFilt(X,N1)

% Xsm = rectFilt(X,N1,N2)
% 
% 1-d rectangular filtering with border correction
% Input:
%     X: vector (or column matrix) to smooth
% Output:
%     Xsm: smoothed matrix

% Adrien Peyrache 2011

if length(size(X))>2
    error('Too many dimensions!')
end

l1 = size(X,1);
X = [flipud(X);X;flipud(X)];

gw = ones(N1,1); %alpha = 0.05
gw = gw/sum(gw);
X = convn(X,gw,'same');
X = X(l1+1:2*l1,:);

