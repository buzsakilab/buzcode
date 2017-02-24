function H = toHsv(data,occ,varargin)

% H = toHsv(data,occ)
% transforms matrix data into HSV image where brightness codes for values of the occ matrix
% data and occ are required to be of the same size
% H = toHsv(data,occ,[minVal maxVal])
% will normalize the values of 'data' matrix in the boundaries define by [minVal maxVal]
% 
% Adrien Peyrache, 2012


H = ones(size(data,1),size(data,2),3);

if isempty(varargin)
    H(:,:,1) = 0.6667*(1-contrast2(data));
else
    rg = varargin{1};
    H(:,:,1) = 0.6667*(1-contrast2(data,rg(2),rg(1)));
end
%H(:,:,3) = 0.9*contrast(occ)+0.1;
H(:,:,3) = 0.9*contrast2(occ)+0.1;
H(:,:,2) = 1;
