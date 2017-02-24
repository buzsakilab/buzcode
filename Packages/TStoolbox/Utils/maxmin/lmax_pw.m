function [lmval,indd]=lmax_pw(xx, dx)
%   Find  piece-wise  local maxima in vector XX,where
%	LMVAL is the output vector with maxima values, INDD  is the 
%	corresponding indexes, DX is length of piece where maxima is searched, 
%   IMPORTANT:      FIRST and LAST point in vector are excluded
%   IMPORTANT:      XX must be single column vector
%   IMPORTANT:      Length of DX must be very carefully selected 
%	For example compare dx=30; and dx=1000;
%
%   dx=150; xx=[0:0.01:35]'; y=sin(xx .* cos(xx /4.5)) + cos(xx); 
%    plot(xx,y); grid; hold on;
%   %   Excluding first and last points
%   [b,a]=lmax_pw(y,dx); plot(xx(a),y(a),'r+')
%   % Way to include first and last points can be as:
%   y(1)=1.5; yy=[0; y; -1;];   % padd with smaller values
%   [b,a]=lmax_pw(yy,dx); a=a-1; plot(xx(a),y(a),'go')
%
%	see also LMIN, LMAX,  LMIN_PW, MATCH

% 	Sergei Koptenko, Applied Acoustic Technologies, Toronto, Canada
%   sergei.koptenko@sympatico.ca,  March/11/2003  

if nargin <2, 
	disp('Not enough arguments'); return
end

len_x = length(xx);
xx = [xx; xx(len_x); xx(len_x)]; 
nn=floor(len_x/dx);
ncount=1; lmval=[]; indd=[];
	for ii=1:nn,
        [lm,ind] = max(xx(ncount: ii*dx+2)) ;
        ind=ind+(ii-1)*dx;
                 if (ind ~=ncount) & (ind~=ii*dx+2),    
                    lmval=[lmval, lm]; indd=[indd, ind]; 
                end      
        ncount=ncount +dx;
	end
[lm,ind] = max(xx(ii*dx:len_x));
        if (ind ~=len_x) & (ind~=ii*dx),    
            lmval=[lmval, lm]; indd=[indd, (ind+ii*dx-1)]; 
        end
    
       if indd(end)==len_x,  
           indd=  indd(1:end-1); 
           lmval=lmval(1:end-1);    
       end
return
