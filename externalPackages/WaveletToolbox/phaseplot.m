function h = PhasePlot(x,y,phase,sz,varargin)
% Plots arrows indicating the phase on current axis.
%
% [h,yy,zz] = PhasePlot(x,y,phase,sz[,ArrowParameters])
% 
% 
% Requires the Arrow function by Erik A. Johnson.
%
% Note: Arrows will get skewed if you resize the axis.
%
% Aslak Grinsted 2003


% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


x=x(:);
y=y(:);
phase=phase(:);
if (length(x)*length(y)==length(phase)) 
    [x,y]=meshgrid(x,y);
end;
sz=sz(:);

%remove nans
idx=find(~any(isnan([x(:) y(:) phase]),2));
x=x(idx);
y=y(idx);
phase=phase(idx);
if (length(sz)>1), sz=sz(idx); end;


dar=get(gca,'DataAspectRatio');
%len=sqrt(dar(1).^2+dar(2).^2);
dxlim=abs(diff(get(gca,'xlim')));
dar(2)=dar(2)*dxlim/dar(1);
dar(1)=1*dxlim;


dx=cos(phase).*dar(1).*sz*.5;
dy=sin(phase).*dar(2).*sz*.5;

h=arrow([x-dx y-dy],[x+dx y+dy],varargin{:});
set(h,'clipping','on')
if (nargout<1) clear h; end;


