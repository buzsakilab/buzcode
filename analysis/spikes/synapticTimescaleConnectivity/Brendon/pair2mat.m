function M = pair2mat(data,id,varargin);

% pair2mat - transform a serie of (x,y) value and associated data to a matrix
% 
% USAGE
%     M = pair2mat(data,id);
%     
%     data: a vector (or 2 column matrix) of values to fill in the output
%           matrix. If data is a vector, the matrix will be symetric. If
%           not, the first (second) column will fill the upper (lower)
%           triangle of the matrix.
%     id: a 2 column matrix of element indices. Must be the same length as data
%
%     When no other argument is secified, the output matrix has the same
%     size as the maximal index of 'id'
% 
%     M = pair2mat(data,id,dim);
%     same usage, but the outt matrix has now a size [dim,dim]
%     dim must be greater than the maximal index in 'id'
% 
%     M = pair2mat(data,id,M0);
%     the matrix 'M0' will be updated with the values in 'data'
%     

% Adrien Peyrache, 2014


dim=0;
if isempty(varargin)
    dim = max(id(:));
else
    if max(size(varargin{1})) == 1
        dim = varargin{1};
    else
        data0 = varargin{1};
        if size(data0,1) ~= size(data0,2)
            error('Original matrix is not squared')
        elseif size(data0,1) < max(id(:))
            error('Original matrix is smaller than some indices')
        end
    end
end

xx = id(:,1);
yy = id(:,2);

if dim==0
    M = data0;
else
    M = NaN(dim,dim);
end
for ii=1:size(xx,1)
    M(xx(ii),yy(ii)) = data(ii,1);
    M(yy(ii),xx(ii)) = data(ii,2);
end
