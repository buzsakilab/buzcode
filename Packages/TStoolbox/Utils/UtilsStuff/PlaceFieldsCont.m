function [pf,x1,x2] = PlaceFieldsCont(S,X,Y,epoch,varargin);

% Place-Field for continuous values
% USAGE
%  [pf,xb,xb] = PlaceFieldCont(S,X,Y,epoch,...);
%  
%  Input:
%  	S: a tsd
%  	X,Y: tsd od spatial position
%  	epoch: intervalset on which place field is computed
%  
%  	by default, the number of bins in the 2 dimensions is 100. You can specify it in two ways
%  	... = PlaceFieldCont(S,X,Y,epoch,bins);
%  	where bins is the number of bins (the same for the 2 dimensions)
%  	... = PlaceFieldCont(S,X,Y,epoch,Xbins,Ybins);
%  	for different number of bins along X and Y
%  
%  Output:
%  	pf: placefield matrix
%  	xb: number of bins along X
%  	yb: number of bins along Y
%
% Adrien Peyrache, 2008

XS = Restrict(X, epoch);
YS = Restrict(Y, epoch);
S = Restrict(S,epoch);

l=length(varargin);

if l==0
	
	bins = 100;
	x1 = bins;
	x2 = bins;

elseif l==1

	bins = varargin{1};
	if ~(max(size(bins))==1 & isnumeric(bins))
		error('Bins should be a scalar')
	end
	x1 = bins;
	x2 = bins;

elseif l==2

	x1 = varargin{1};
	x2 = varargin{2};
	if ~(max(size(x1))==max(size(x2)) & isnumeric(x1) & isnumeric(x2))
		error('Problem with X1 and X2')
	else
		bins = length(x1);
	end
else
	error('Problem with your arguments!!')
end

    %db = round(bins/10);

    pX = Restrict(XS, S);
    pY = Restrict(YS, S);
    pf = histCont2d(Data(S),Data(pX), Data(pY), x1, x2);

    %largerMatrix = zeros(bins+2*db,bins+2*db);
    %largerMatrix(db+1:bins+db,db+1:bins+db) = pf;
    %pf = largerMatrix;

    pf(isnan(pf))=0;
    %gw = gausswin(round(db));
    %gw = gw*gw';
    %gw = gw/sum(sum(gw));

    %pf = convn(pf,gw,'same');

