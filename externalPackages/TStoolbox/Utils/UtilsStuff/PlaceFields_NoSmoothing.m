function [pf,occH,bx,by,pX,pY] = PlaceField(S,X,Y,epoch,varargin)

% [pf,occH,x1,x2] = PLACEFIELD(tsa,X,Y,epoch,bins,occH,x1,x2);
%  Different Usages
%  	1- PlaceField(tsa,X,Y,epoch,bins);
%  		where bins is the size of the bins in cms
%  	2- PlaceField(tsa,X,Y,epoch,x1,x2);
%  		you specify the binning spacing x1,x2 along X and Y respectively
%       OR the number of bins along each dimensions
%  	3- PlaceField(tsa,X,Y,epoch,x1,x2,occH);
%  		you have already computed the occupancy matrix occH, no need to loose time,
%  		give it to the program.
%  OUTPUT:
%  		pf: the placefield map
%  		occH: occupancy map
%  		x1: bins spacing along the X
%  		x1: bins spacing along the Y
%  
%  Adrien Peyrache 2008


XS = Restrict(X, epoch);
YS = Restrict(Y, epoch);
S = Restrict(S,epoch);

l=length(varargin);

if l==1

	bins = varargin{1};
	if ~(max(size(bins))==1 & isnumeric(bins))
		error('Bins should be a scalar')
    else
        x = Data(XS);
        y = Data(YS);
        nbBinsX = round((max(x)-min(x))/bins);
        nbBinsY = round((max(y)-min(y))/bins);

        %Inverted to preserve (x,y)
		[occH, by, bx] = hist2([y, x], nbBinsX, nbBinsY);
    end

elseif l==2

	bx = varargin{1};
	by = varargin{2};
    [occH by bx] = hist2([Data(YS), Data(XS)], by, bx);
    nbBinsX = length(bx);
    nbBinsY = length(by);
    
elseif l==3
	bx = varargin{1};
	by = varargin{2};
	occH = varargin{3};

	if ~(max(size(bx))>1) | ~max(size(by))>1
		error('Problem with X1 and X2')

	elseif ~(size(occH,1)==max(size(by)) & size(occH,2)==max(size(bx)) & isnumeric(occH))
		error('Problem with X1 and X2')
	else
        nbBinsX = length(bx);
        nbBinsY = length(by);
	end
else
	error('Problem with your arguments!!')
end

pX = Restrict(XS, S);
pY = Restrict(YS, S);

if length(pX)

    pfH = hist2([Data(pY), Data(pX)], by, bx);

    warning off
    pf = 39.0625 * pfH./occH;
    warning on

    pf(isnan(pf))=0;
    pf(isinf(pf))=0;

else
    pf = zeros(nbBinsY,nbBinsX);
end
