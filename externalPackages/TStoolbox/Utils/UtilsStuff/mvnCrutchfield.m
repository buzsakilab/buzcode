function d = mvnCrutchfield(Data1,Data2,varargin)

%  multivariate crutchfield distance computation
%  Adrien Peyrache 2007


N = size(Data1,2);
if size(Data2,2) ~= N
	error('Data dimensions don''t fit')
end

if nargin ==1
	mask = varargin{1};
	if (size(mask,1)~=N) && (size(mask,2) ~=N)
		error('mask dimensions don''t fit')
	end
elseif nargin>1
	warning('too many input arguments');
end


nbS1 = size(Data1,1);
nbS2 = size(Data2,1);

Data = [Data1,Data2];
c = corrcoef(Data);
c11 = c(1:N,1:N);
c22 = c(N+1:end,N+1:end);
c12 = c(N+1:end,1:N);
c21 = c(1:N,N+1:end);

cCond = c11 - c12*inv(c22)*c21;
H1 = N + log(det(cCond));

Data = [Data2,Data1];
c = corrcoef(Data);
c11 = c(1:N,1:N);
c22 = c(N+1:end,N+1:end);
c12 = c(N+1:end,1:N);
c21 = c(1:N,N+1:end);

cCond = c11 - c12*inv(c22)*c21;
H2 = N + log(det(cCond));

keyboard

d = H1+H2;