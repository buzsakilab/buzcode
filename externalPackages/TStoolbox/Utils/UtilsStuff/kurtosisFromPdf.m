function k = kurtosisFromPdf(X,P)

%  k = kurtosisFromPdf computes the kurtosis of the distribution given by the pdf P
%  over the vector of data X
%  
%  example: 
%  X = [-5:5];
%  pdf = 1/sqrt(2*pi) * exp(-(X.^2)/2);
%  k = kurtosisFromPdf(X,pdf)
%  
%  gives:
%  
%  k =
%  
%      3.0000
%  
% as expected for a normal distribution
% Adrien Peyrache 2007


a = size(X);
[n,ix] = max(a);



if ix > 2
	error('X must be a vector')
elseif ix==1
	if size(X,2)>1
		error('X must be a vector')
	end
	
	X = X';
	P = P';
end;

n = length(X);
P=P./(sum(P));

mu = X*P';
x2 = ((X-mu).*(X-mu))*P';
x4 = ((X-mu).*(X-mu).*(X-mu).*(X-mu))*P';

k = x4/(x2^2);
%  keyboard;


%  keyboard


