function [lmval,indd]=lmin_pw(xx, dx)
%   Find  piece-wise local minima in vector XX,where
%	LMVAL is the output vector with minima values, INDD  is the 
%	corresponding indexes, DX is scalar length of piece where minima is searched, 
%   IMPORTANT:      FIRST and LAST point in vector are excluded
%   IMPORTANT:      XX must be single column vector
%   IMPORTANT:      Length of DX must be very carefully selected 
%	For example compare dx=10; and dx=1000;
%
%   dx=150; xx=[0:0.01:35]'; y=sin(xx .* cos(xx /4.5)) + cos(xx); 
%   y(length(y))=-2; plot(xx,y); grid; hold on;
%   %   Excluding first and last points
%   [b,a]=lmin_pw(y,dx); plot(xx(a),y(a),'r+')
%   % Way to include first and last points can be as:
%   yy=[1.5; y; 0];         % padd with values higher than end values
%   [b,a]=lmin_pw(yy,dx); a=a-1; plot(xx(a),y(a),'go')
%
%	see also LMIN,LMAX, LMAX_PW, MATCH

% 	Sergei Koptenko, Applied Acoustic Technologies, Toronto, Canada
%   sergei.koptenko@sympatico.ca,  March/11/2003  

if nargin <2, 
	disp('Not enough arguments'); return
end

len_x = length(xx);
xx = [xx; xx(end)]; 
nn=floor(len_x/dx);
ncount=1; lmval=[]; indd=[];
for ii=1:nn,
    [lm,ind] = min(xx(ncount: ii*dx+1)) ;
        ind=ind+(ii-1)*dx;
         if (ind ~=ncount) & (ind~=ii*dx+1),    
         lmval=[lmval, lm]; indd=[indd, ind]; 
end      
ncount=ncount +dx;
end
[lm,ind] = min(xx(ii*dx:len_x));
    if (ind ~=len_x) & (ind~=ii*dx),    
    lmval=[lmval, lm]; indd=[indd, (ind+ii*dx-1)]; 
    end
    
     if indd(end)==len_x,
    indd=indd(1:end-1); lmval=lmval(1:end-1);
    end
return
