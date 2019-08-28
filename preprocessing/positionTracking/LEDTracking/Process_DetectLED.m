% Detect red and blue LEDs position in a video file and creates a 'led' file
%
%  USAGE
%
%    Process_DetectLED(videoFile,<options>)
%
%    videoFile      path to Basler video file, including the '.avi'
%                   extension or not
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'manualROI'   boolean: whether or not you want to manually adjust the
%                   area where LEDs are detected (default = 1)
%     'ROI'         2 X C matrix where C1 is x and C2 is y from ROI
%    =========================================================================
%
% DEPENDENCIES:
%
%   Computer vision toolbox


% Copyright (C) 2015 Adrien Peyrache, some inspiration from John D Long II
% 2019 Manu Valero (added roi inputs)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


function whl = Process_DetectLED(fbasename,thresh,varargin)
% Parse options
p = inputParser;
addParameter(p,'manualROI',false,@islogical);
addParameter(p,'roi',[],@ismatrix);
parse(p,varargin{:});

manualROI = p.Results.manualROI;
roi = p.Results.roi;

if strcmp(fbasename(end-3:end),'avi')
    fbasename = fbasename(end-3:end);
end
% 
file = [fbasename '.avi'];
% file = fbasename;

if ~exist(file,'file')
    warning('No video file')
    keyboard
end

% Creater readerobj of file
videoObj    = vision.VideoFileReader(file);
videoSize   = videoObj.info.VideoSize;
width       = videoSize(1);
height      = videoSize(2);
threshF     = thresh;% 0.15;    % threshold on foreground pixel intensity

% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);

% Define ROI (manual definition of the environment)
%try 
%First, we adapt the dynamical range of the pixel colors for manual
%selection of the ROI (i.e. the environment)
firstFrame  = step(videoObj);
frame = imadjust(rgb2gray(firstFrame));

if manualROI
    ok = 0;
    figure(1),clf

    while ~ok
        clf,imshow(frame)
        fprintf('Define your ROI. Click ''enter'' when finished\n')
        [x,y] = ginput;
        inArea = inpolygon(X(:),Y(:),x,y);
        inArea = double(reshape(inArea,[height width]));
        frame(~inArea) = 0;

        clf,imshow(frame)

        reply = input('OK with the result? Y/N [Y]:','s');
        if ~strcmp(reply,'N') || ~strcmp(reply,'n')
            ok = 1;
        end
    end
elseif ~manualROI && ~isempty(roi)
    x = roi(:,1);
    y = roi(:,2);
    inArea = inpolygon(X(:),Y(:),x,y);
    inArea = double(reshape(inArea,[height width]));
    frame(~inArea) = 0;
else
    inArea = ones(size(X));
end
% Initialize background as a grayscale image of the first frame
bg_bw       = rgb2gray(firstFrame);

% Initialize foreground image and difference image
fg          = zeros(size(bg_bw));
fr_diff     = zeros(size(bg_bw));

% Initialize color mask
mask  = zeros(height,width,3,'uint8');

% Initialize fr in case first frame reading returns error
fr  = zeros(height,width,3,'uint8');
bkgound = firstFrame;
% Initialize whl matrix
whl = [];
count = 0;

try   
while ~isDone(videoObj)
%for i = Fint:readerobj.NumberOfFrames
    
    if count~=0
        backSp = repmat('\b',[1 length(num2str(count-1))]);
        fprintf(backSp)
    end
    fprintf('%i',count)

    if count~=0
        fr    = step(videoObj);
    else
        fr = firstFrame;
    end
%     bkground = (fr + bkground) ./ 2;
    
    % convert frame to grayscale
    fr_bw = rgb2gray(fr);
    %%% Label Color Mask
    label         = repmat(logical(fr_bw>threshF) & inArea,[1 1 3]);
    mask(label)   = fr(label);
    mask(~label)  = 0;
    
    %%% Find centroid of remaining pixels %%%
    %Red
    bw_mask = squeeze(mask(:,:,1));
    [CC,nc] = bwlabel(bw_mask);
    
    if nc>0
        pixels = regionprops(CC,'PixelList');
        centroidSize = zeros(nc,1);
        for ii=1:nc
            p = pixels(ii).PixelList;
            centroidSize(ii) = length(p);
        end
        [~,mxIx] = max(centroidSize);
        Rr   = round(mean(pixels(mxIx).PixelList,1));
    else
        Rr = [-1 -1];
    end
    
    %Blue
    Br = [-1 -1];
    bw_mask = squeeze(mask(:,:,3));
    [CC,nc] = bwlabel(bw_mask);
    pixels = regionprops(CC,'PixelList');
    
    if nc>0
        pixels = regionprops(CC,'PixelList');
        centroidSize = zeros(nc,1);
        for ii=1:nc
            p = pixels(ii).PixelList;
            centroidSize(ii) = length(p);
        end
        [~,mxIx] = max(centroidSize);
        Br   = round(mean(pixels(mxIx).PixelList,1));
    else
         Br = [-1 -1];
    end
    
     whl = [whl;[Rr(1),Rr(2),Br(1),Br(2)]];
    
    % End processing time
    
    % Display results every 1000 frame 
    if 0 %set it at 1 for debug
        if mod(count,1000)==0
%             figure(1),h1 = subplot(3,1,1);imshow(uint8(fr)), %title(sprintf('Frame %d, Time %f, Connected Components %d',i,t2,Nl))
%             subplot(3,1,2),imshow(uint8(floor(fr_diff)))
%             h2 = subplot(3,1,3);imshow(mask)
            %h3 = impoly(h1,[xj,yj]);
            %h4 = impoly(h2,[xi,yi]);
            %setColor(h3,'green')
            %setColor(h4,'yellow')
            %keyboard
            figure(2),
            plot(whl(:,1),whl(:,2),'r')
            hold on
            plot(whl(:,3),whl(:,4),'b')
            drawnow
        end
    end
    count = count+1;
end
catch
   keyboard
end
%write result file
dlmwrite([fbasename '.led'],whl,'\t')

fprintf('\n\n')
release(videoObj)