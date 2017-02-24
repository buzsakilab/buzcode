function [lmval,indd]=lmin(xx,filt)
%LMIN 	function [lmval,indd]=lmin(x,filt)
%	Find local minima in vector X, where LMVAL is the output
%	vector with minima values, INDD is the corresponding indeces 
%	FILT is the number of passes of the small running average filter
%	in order to get rid of small peaks.  Default value FILT =0 (no
%	filtering). FILT in the range from 1 to 3 is usially sufficient to 
%	remove most of a small peaks
%	Examples:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
%	plot(xx,y); grid; hold on;
%	[a b]=lmin(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMAX, MAX, MIN
	
%
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
	fltr=[1 1 1]/3;
  if nargin <2, filt=0; 
	else
x1=x(1); x2=x(len_x); 

	for jj=1:filt,
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end
  end

lmval=[];
indd=[];
i=2;		% start at second data point in time series

    while i < len_x-1,
	if x(i) < x(i-1)
	   if x(i) < x(i+1)	% definite min
lmval =[lmval x(i)];
indd = [ indd i];

	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case 
%indd = [ indd i];	%2 when only  definite min included
i = i + 2;  		% skip 2 points

	   elseif x(i)==x(i+1)	% 'short' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite min included
i = i + 1;		% skip one point
	   end
	end
	i = i + 1;
    end

if filt>0 & ~isempty(indd),
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
	   rng=1;	%check if index too close to the edge
	else rng=2;
	end

	   for ii=1:length(indd), 
		[val(ii) iind(ii)] = min(xx(indd(ii) -rng:indd(ii) +rng));
		iind(ii)=indd(ii) + iind(ii)  -rng-1;
	   end
  indd=iind; lmval=val;
else
end

