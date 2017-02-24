function [pf,occH,x1,x2,pX,pY] = PlaceField(S,X,Y,epoch,varargin);

% [pf,occH,x1,x2] = PLACEFIELD(tsa,X,Y,epoch,bins,occH,x1,x2);
%  Different Usages
%  	1- PlaceField(tsa,X,Y,epoch,bins);
%  		where bins is the number of bins on X and Y
%  	2- PlaceField(tsa,X,Y,epoch,x1,x2);
%  		you specify the binning spacing x1,x2 along X and Y respectively
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
smoothOcc = 1;

l=length(varargin);

if l==1

	bins = varargin{1};
	if ~(max(size(bins))==1 & isnumeric(bins))
		error('Bins should be a scalar')
    else
        display(bins)
        %Inverted to preserve x,y
		[occH, x2, x1] = hist2([Data(YS), Data(XS)], bins, bins);
        smoothOcc = 1;
    end

elseif l==2
	x1 = varargin{1};
	x2 = varargin{2};
	if ~(max(size(x1))==max(size(x2)) & isnumeric(x1) & isnumeric(x2))
		error('Problem with X1 and X2')
	else
		[occH, x2, x1] = hist2([Data(YS), Data(XS)], x2, x1);
		bins = length(x1);
        smoothOcc = 1;
	end

elseif l==3
	x1 = varargin{1};
	x2 = varargin{2};
	occH = varargin{3};

	if ~(max(size(x1))==max(size(x2)) & isnumeric(x1) & isnumeric(x2))
		error('Problem with X1 and X2')

	elseif ~(size(occH,1)==max(size(x2)) & size(occH,2)==max(size(x2)) & isnumeric(occH))
		error('Problem with X1 and X2')
	else
		bins = length(x1);
	end
else
	error('Problem with your arguments!!')
end

db = round(bins/5);

    pX = Restrict(XS, S);
    pY = Restrict(YS, S);

    if length(pX)
	pfH = hist2([Data(pY), Data(pX)], x2, x1);

	
% 	keyboard
	%pfH = gausssmooth(pfH,db,db);   
    %occH = gausssmooth(occH,db,db); 
    dt = 1/median(diff(Range(X,'s')));
    

    warning off
	pf = dt * pfH./occH;
    warning on

	%largerMatrix = zeros(bins+2*db,bins+2*db);
	%largerMatrix(db+1:bins+db,db+1:bins+db) = pf;
	%pf = largerMatrix; % X and Y axis are inverted, so need to take the transpose;
	
	pf(isnan(pf))=0;
	pf(isinf(pf))=0;
	pf = gaussFilt(pf,db/5);   
    else
        pf = zeros(bins+2*db,bins+2*db);
    end

    if 0PlaceFields
        
        figure(1),clf
        subplot(1,2,1)
        imagesc(x2,x1,pf),axis xy
        subplot(1,2,2)
        plot(Data(XS),Data(YS),'Color',[0.5 0.5 0.5])
        hold on
        plot(Data(pX),Data(pY),'r.')
        
    end