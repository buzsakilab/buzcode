function [Xout,Yout]=gplot3(A,xy,lc)

[i,j] = find(A);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-separated list of line segments,
% rather than individual segments.

X = [ xy(i,1) xy(j,1) repmat(NaN,size(i))]';
Y = [ xy(i,2) xy(j,2) repmat(NaN,size(i))]';
Z = [ xy(i,3) xy(j,3) repmat(NaN,size(i))]';
X = X(:);
Y = Y(:);
Z = Z(:);

if nargout==0,
    if nargin<3,
        plot3(X, Y, Z)
    else
        plot3(X, Y, Z, lc);
    end
else
    Xout = X;
    Yout = Y;
end