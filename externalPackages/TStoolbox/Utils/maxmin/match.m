function indd=match(x,y, tol,bro)
%   MATCH  function INDD=MATCH(X, Y, TOL, BRO) find the closest
%	values to X in vector Y with tolerance TOL. If TOL  is omitted, the
%   default  value of 1e-6 will be used. CAUTION: Using too small tolerance may
%   exclude some points from results (run example below with tol=0.005)
%   BRO=0 by default and can also be omitted
%   If BRO=0, then all indices of values of Y within the tolerance TOL vill be
%   included in the output. If BRO=1 then one (last or first) value within each range
%   will be reported.  If BRO>1 then one minimum value for each range
%   will be included in output. Value of BRO defines indices range for
%   minimum matching value search (rule of thumb- the coarser your TOL 
%   the bigger BRO should be. Use BRO= twice the maximum number of
%   matching values for each matching point)
%   (Run example below with a=match(val,y,tol, 1);
%   and tol=0.01; to compare outputs)
%	For example:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
%	plot(xx,y); grid; hold on;
%	val=1.04; tol=0.05;
%	a=match(val,y,tol); plot(xx(a),y(a),'r.');  % All matching values
%   aa=match(val,y,tol,1);  plot(xx(aa),y(aa),'ms'); % Single value per range
%   aaa=match(val,y,tol,20);  plot(xx(aaa),y(aaa),'go'); % Min value
%
%	see also LMIN,LMAX, LMAX_PW,  LMIN_PW

% 	Sergei Koptenko, Applied Acoustic Technologies, Toronto, Canada
%   sergei.koptenko@sympatico.ca,  March/22/2003  

	if (nargin<4),       bro=0; 
	else,    bro=round(bro);
	end;		% set the default  broadness
		if (nargin<3),   tol=1e-6;  bro=0;  % set the default  tolerance
        end;		
        
 [rrow,ccol]=size(y);
if  (ccol==1), y=y'; end
[rrow,ccol]=size(y);
if  (rrow>1) || (rrow+ccol==2), disp('ERR: Y must be a vector!'); indd=[]; return; end
if ( x>max(y)) || (x< min(y) ), 
    disp('ERR: X is outside of range of Y!');
    indd=[]; return; 
end

indd= find(xor(((x+tol) >=y), ((x-tol)<=y)) ==0); %get all indices within the tolerance

if bro>0,
	difi= diff(indd);              % find index difference
	fdif= find(~(difi-1));     % find neigbouring indices 
	ddff=indd(fdif);           % get neigbouring indices
%----------remove neigbouring indices-------------------
				for k1=1:length(ddff),
                    aa = find(~(indd - ddff(k1)));
                    indd(aa)=0;
				end %k1
	bb=sort(indd);
	indd=bb(find(bb));
%----------find closest value for each range-------------------
			if bro>1,
                innd=1;
                        for j3=1:length(indd),
                               j4=indd(j3)-bro;
                               j5=indd(j3)+bro;
                               dj4=0;
                                    if j4<1, dj4=-j4 +1;
                                        j4=1;
                                    end
                                    if j5>ccol, j5=ccol;end
                               [vv,innd(j3)] = min(abs( y(j4:j5) -x ));
                               innd(j3)=innd(j3) + dj4;
                         end %j3          
                            indd=indd + innd - bro-1;
             end
end 
                          if ~isempty(indd), 
                          if indd(1)<1, indd(1)=1;  end
                           if indd(end)>ccol, indd(end)=ccol; end
                       end
return
