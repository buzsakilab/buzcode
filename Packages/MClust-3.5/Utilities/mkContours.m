function mkContours(axesFrom, varargin)
% make contour plot
% 
% Original code by Peter Lipa 1999
%
% Modified by ADR 2000
%
% INPUTS
%   axesFrom = axis handle from which to draw the scatter points
% 
% parameters
%   nX = 100
%   nY = 100
%   figHandle = figure to draw in (def = new fig)
%   
%
% Draws a contour plot from scatter plot
%

% ADR 2000
%
% Status: PROMOTED (Release version) 
% See documentation for copyright (owned by original authors) and warranties (none!).
% This code released as part of MClust 3.0.
% Version control M3.0.
% adr 17 aug 2003 added toolbox warning

nX = 100;
nY = 100;
figHandle = [];
extract_varargin;

% -- get axis
fhaxes = findobj(axesFrom, 'Type', 'Axes');

%--- get some axis properties
XLim = get(fhaxes,'XLim');
YLim = get(fhaxes,'YLim');
AxLabels{1} = get(get(fhaxes,'XLabel'),'String');
AxLabels{2} = get(get(fhaxes,'YLabel'),'String');

%--- find figure handles to cluster line objects (identified by linestyle = none)
fhc = findobj(fhaxes, 'Type', 'Line'); 
%--- get and merge XData and YData
XData = [];
YData = [];
for ic = 1:length(fhc);
   Xc = get(fhc(ic),'XData');
   Yc = get(fhc(ic),'YData');
   if length(Xc) > 1
      XData = [XData Xc];
      YData = [YData Yc];
   end%if
end%for ic

if isempty(XData) 
  return
end

%--- finally make the histogram
H = ndhist([XData; YData], [nX; nY], [XLim(1); YLim(1)], [XLim(2); YLim(2)])';

S = Hsmooth(H);
lS = fixedLog(S,-1);

if isempty(figHandle)
   figHandle = figure;
end
[nX, nY] = size(lS);
X = linspace(XLim(1),XLim(2),nX);
Y = linspace(YLim(1),YLim(2),nY);
figure(figHandle);
Cont = contourf(X, Y , lS, 15);
xlabel(AxLabels{1}); ylabel(AxLabels{2});
zoom on;
% colorbar;

%==========================================================
function fL = fixedLog(H, infimum)

%  returns the log(H) with H a N-dim array 
%  where the NaNs of log(x<=0) of H are replaced by infimum
%

fL = H;
fL(fL <0) = 0;
fL = log(fL);
fL(fL <= infimum) = infimum;

return;