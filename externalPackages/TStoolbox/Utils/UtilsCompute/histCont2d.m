function [no,xo1,xo2] = histCont2d(z,y1,y2,x1,x2)
% HISTCONT2D  2-D histogram of continuous measure
%   [N, X] = HIST2(z,Y1,Y2) bins the elements of Y1 and Y2 
%   into 10x10 equally spaced containers
%   and returns the average value of z for each bins.
%
%   N = HIST(z,Y1,Y2,M,N), where M and N are scalars, uses MxN bins.
%
%   N = HIST(z,Y1,Y2,X,Y), where X is a matrix, returns the distribution of Y
%   among bins with centers specified by X.
%
%   [N,X1,X2] = HIST(...) also returns the position of the bin centers in X1,X2.
%
%   HIST(...) without output arguments produces a colormap plot of
%   the results.
%
% adapted by Adrien Peyrache, 2008



if nargin < 2
    error('Requires  at least two  input arguments.')
end
if nargin == 3
    x1 = 10;
    x2 = 10;
end
if min(size(y1))==1, y1 = y1(:); end
if min(size(y2))==1, y2 = y2(:); end
% if isstr(x1) | isstr(x2) | isstr(y1) | isstr(y2)
%     error('Input arguments must be numeric.')
% end

[m,n] = size(y1);
if ~all(size(y2) == [m,n])
    error('Y1 and Y2 must have the same size.')
end
if length(x1) == 1
    miny1 = min(min(y1));
    maxy1 = max(max(y1));
    binwidth1 = (maxy1 - miny1) ./ x1;
    xx1 = miny1 + binwidth1*(0:x1);
    xx1(length(xx1)) = maxy1;
    x1 = xx1(1:length(xx1)-1) + binwidth1/2;
else
    xx1 = x1(:)';
    miny1 = min(min(y1));
    maxy1 = max(max(y1));
    binwidth1 = [diff(xx1) 0];
    xx1 = [xx1(1)-binwidth1(1)/2 xx1+binwidth1/2];
    xx1(1) = miny1;
    xx1(length(xx1)) = maxy1;
end

if length(x2) == 1
    miny2 = min(min(y2));
    maxy2 = max(max(y2));
    binwidth2 = (maxy2 - miny2) ./ x2;
    xx2 = miny2 + binwidth2*(0:x2);
    xx2(length(xx2)) = maxy2;
    x2 = xx2(1:length(xx2)-1) + binwidth2/2;
else
    xx2 = x2(:)';
    miny2 = min(min(y2));
    maxy2 = max(max(y2));
    binwidth2 = [diff(xx2) 0];
    xx2 = [xx2(1)-binwidth2(1)/2 xx2+binwidth2/2];
    xx2(1) = miny2;
    xx2(length(xx2)) = maxy2;
end
nbin1 = length(xx1);
nbin2 = length(xx2);
%  xx1 = [0 xx1];
%  xx2 = [0 xx2];
nn = zeros(nbin1-1,nbin2-1,n);

%  keyboard

for i=2:nbin1
   for j=2:nbin2
      t = y1 < xx1(i) & y1 >= xx1(i-1) & y2 < xx2(j) & y2 >=  xx2(j-1); 
      if sum(t)
        nn(i-1,j-1,:) = mean(z(t));
      else
        nn(i-1,j-1,:) = 0;
      end
  end
end

    no = nn';
    xo1 = x1;
    xo2 = x2;

