function [Idx,C] = ClusterSpkTrain(data,n)

%  [Idx,C] = ClusterSpkTrain(data,n) clusters binned spike trains in data into n clusters
%  
%  INPUT : 
%  	data: a matrix in the form time bins x cells
%  	n: number of clusters
%  OUTOUT:
%  	Idx: cluster indices associated with each cell
%  	C: simimlarity matrix
%
%  adapted from J-M Fellous, J Neuroscience 2004
%  Adrien Peyrache 2007,2009


l= size(data,2);
C = eye(l);

C = data'*data;

normD = sqrt(sum(data.^2));
C = C./(normD'*normD);

diagIx = eye(size(C));
c = C(:);
s = c(~diagIx);
meanS = mean(mean(s));

tauVal = [0.01:0.01:0.3];
sdB = [];

for tau=tauVal
	
	b = 1./(1+exp((meanS-s(:))/tau));
	[H,B] = hist(b);
	sdB = [sdB;std(H)];

end

[m,ix] = min(sdB);
tau = tauVal(ix);

b = 1./(1+exp((meanS-c(:))/tau));
c = reshape(b,size(C));

[Idx,Centers] = kmeans(c,n,'distance','correlation');

[val,ix] = sort(Idx);
cM = c-diag(diag(c));

if 0

	figure(1),clf
	imagesc(C)

	figure(2),clf
	imagesc(cM)

	figure(3),clf
	imagesc(cM(ix,ix))


end
