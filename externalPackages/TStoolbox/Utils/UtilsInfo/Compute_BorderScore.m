% Compute border score as in Solstad 2008 (works with square matrix)
%
%  USAGE
%
%    Compute_BorderScore(pf,bx,by,<options>)
%
%    pf             place field (#by x #bx)
%    bx             x-axis bins
%    by             y-axis bins
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'manualROI'   boolean: whether or not you want to manually adjust the
%                   area where LEDs are detected (default = 1)
%    =========================================================================
%
% DEPENDENCIES:
%
%   Computer vision toolbox


% Copyright (C) 2015 Adrien Peyrache, some inspiration from John D Long II
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

function [borderSc wallBorder] = Compute_BorderScore(pf,bx,by,varargin)

minSurf = 200;
percPkRate = 0.5;

if ~isempty(varargin)
    minSurf = varargin{1};
    if length(varargin)>1
        percPkRate = varargin{2};
    end
end
% 
% if size(pf,1) ~= size(pf,2)
%     error('Matrix must be square')
% end

[xm,ym] = meshgrid(bx,by);
areaBin = median(diff(bx))*median(diff(by));

distM = zeros([size(pf),4]);
distM(:,:,4) = ym-by(1); %Distance to south wall
distM(:,:,2) = sqrt((by(end)-ym).^2); %Distance to north wall
distM(:,:,3) = xm-bx(1); %Distance to west wall
distM(:,:,1) = sqrt((bx(end)-xm).^2); %Distance to east wall
dist = squeeze(min(distM,[],3));
wallIx = zeros(size(pf));

wallIx(:,end) = 1; %East
wallIx(end,:) = 2; %North
wallIx(:,1) = 3; %West
wallIx(1,:) = 4; %South

%distSort = sort(dist(:),'descend');

pf = pf/max(pf(:));
[CC,nc] = bwlabel(double(pf>percPkRate));

maxDist = min(bx(end)-bx(1),by(end)-by(1))/2;

borderSc = [];
borderScCt = [];
area = [];
pkRate = [];
wallBorder = [];

if nc>0
    
    fields = regionprops(CC,'PixelList','area');
    
    for ii=1:nc
        a = fields(ii).Area;
        area = [area;a];
        
        wIx = zeros(4,1);
        for w=1:4
            ix = CC==ii & wallIx==w;
            wIx(w) = sum(ix(:))/sum(sum(wallIx==w));
        end
        
        if a*areaBin > minSurf & any(wIx)>0
            p = fields(ii).PixelList;
            I = sub2ind(size(dist),p(:,2),p(:,1));

            [cM,ix] = max(wIx);
            wallBorder = [wallBorder;ix];
            d = squeeze(distM(:,:,ix));
            d = d(I);
            pfI = pf(I);
            pfI = pfI/sum(pfI(:));
            d = pfI'*d/maxDist;
            
            b = (cM-d)/(cM+d);
            borderSc = [borderSc;b];
            pkRate = [pkRate;max(pf(I))];
            
        else
            pkRate = [pkRate;NaN];
            borderSc = [borderSc;NaN];
            wallBorder = [wallBorder;NaN];
        end
        
    end
 
    [~,ix] = max(borderSc);
    borderSc = borderSc(ix);
    area = area(ix);
    wallBorder = wallBorder(ix);
    
    if 0 %borderSc>0.4
        disp([area*areaBin/100 borderSc])
        figure(1),clf
        subplot(1,2,1)
            imagesc(bx,by,pf),axis xy
            colorbar
        subplot(1,2,2)
            imagesc(bx,by,CC==ix),axis xy
        keyboard
    end
else
    borderSc = NaN;

end
