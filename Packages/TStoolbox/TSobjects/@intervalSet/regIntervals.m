function R = regIntervals(ep,N)

%  Split an intervalSet object in N segments of equal length
%  
%  	USAGE
%  	R = regIntervals(ep,N)
%  	
%  	INPUTS:
%  	ep: original intervalSet
%  	N: number of desired segments
%  	
%  	OUTPUT:
%  	R: a cell array of splitted intervalSet objects


% copyright (c) 2004-2011 Francesco P. Battaglia & Adrien Peyrache
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html
 
N =round(N);
R = cell(N,1);

l = Data(length(ep));
ls = sum(l)/N;

st = Start(ep);
en = End(ep);

ix=1;
   
for ii=1:N    
    stR = [];
    enR = []; 
    lR=l(ix);
    while lR<ls & ix<length(l)
        stR = [stR;st(ix)];
        enR = [enR;en(ix)];        
        ix = ix+1;
        lR = lR+l(ix);        
    end
    if lR>ls
        lR = lR-l(ix);
        stR = [stR;st(ix)];
        enR = [enR;st(ix)+ls-lR];
        st = [st(1:ix);st(ix)+ls-lR;st(ix+1:end)];
        en = [en(1:ix-1);st(ix)+ls-lR;en(ix:end)];
        l = [l(1:ix-1);ls-lR;l(ix)-(ls-lR);l(ix+1:end)];
    else
        stR = [stR;st(ix)];
        enR = [enR;en(ix)];
    end
    R{ii} = intervalSet(stR,enR);
    ix=ix+1;
    
end