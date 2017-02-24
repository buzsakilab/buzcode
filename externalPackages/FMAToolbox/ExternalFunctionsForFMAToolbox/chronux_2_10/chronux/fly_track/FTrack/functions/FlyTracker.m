function [x, y, orientation] = FlyTracker(filename, FrameRange,...
                        NBackFrames, opt1, sqrsize,...
                        alpha, ellipse)

%FLYTRACKER
% Usage: 
%  [x, y, orientation] = FlyTracker(filename, FrameRange, ...
%                        NBackFrames, opt1, sqrsize,...
%                        alpha, ellipse)
%
% This function takes an input filename for a movie and tracks the position
% of a fly from StartFrame to EndFrame.  The movie can be in any format
% (.avi, .mpg, .mp2, etc.).  The following are a list of input parameters
% and output variables.
%
% Input:
%      filename:    must be the complete string of the movie filename
%      FrameRange:  range of frames to examine
%      NBackFrames: number of frames in static background, used only to
%                   find fly in first frame
%      opt1:        'neg' or 'pos'.  This multiplies the frame data by a
%                   negative sign ('neg') or not ('pos').  FindFly needs the 
%                   fly to be a bright spot, so the choice of this depends on
%                   the origional image.  If the original image is RGB, chances
%                   are that 'neg' needs to be selected.
%      sqrsize:     This is half the length of the sides of the square 
%                   around the brightest spot to be used in the CM calculation
%                   and the orientation  calculation.
%      alpha:       This is the weighting parameter for the running
%                   background.  It must be close to 1 for best results.
%                   0.9 seems like a good value for this, anything less
%                   seems to degrade SNR.  But for generalities sake, it
%                   is still a user input parameter.
%      ellipse:     parameters of the elliptical area that will be tracked
%
% Output:
%      x,y:         The x and y locations of the fly in each frame      
%      orientation: A matrix containing the UHP and LHP body axis orientation 
%                   angles for each frame.  UHP are in the first row, LHP are in
%                   the second row.  See FlyOrient and FindFly for more
%                   details.

% This function uses the videoReader function to input the video file into a Matlab object.
% Usage on videoReader and other functions in the videoIO toolbox can be
% found in the videoIO documentation.

%Written by Dan Valente 8/8/07

%Please note: Video indexes start at frame 0
StartFrame = FrameRange(1);
EndFrame = FrameRange(end);
NFrames = EndFrame-StartFrame+1;

%If NBackFrames is smaller than the number of frames we want to look at,
% we'll just take a smaller background average using all the frames of
% interest.
if (NBackFrames >= NFrames)
    NBackFrames = NFrames;
end

%initialize some variables
x = [];
y = [];
orientUHP=[];
orientLHP=[];

%load in video object using videoReader
video = videoReader(filename);
info = getinfo(video);
width = info.width;
height = info.height;

%just prevents error if user wants to track to end, but enters the
% incorrect EndFrame (recall frame index starts at 0).
if (info.numFrames == EndFrame)
    EndFrame = EndFrame-1;
elseif (info.numFrames < EndFrame)
    EndFrame = info.numFrames-1;
end

%seek to StartFrame for initial background calculation
seek(video, StartFrame);

mask = ellipse.mask;

disp('Calculating background...')
background = zeros(height, width);
for bframe = 1:NBackFrames; 
    img = getframe(video);
    current = mask.*double(img(:,:,1));  %just take first color matrix and mask out unwanted portions
    background = background + current/NBackFrames ;
    next(video);
end
flipper = 256*ones(height, width).*mask;
if (strcmp(opt1,'neg') == 1)        %make sure our fly will be a bright spot
    background = -background+flipper;
end


%seek to StartFrame for tracking
seek(video,StartFrame);
disp('Beginning tracking...')
tic

for frame = StartFrame:EndFrame
        disp(num2str(frame))
        img = getframe(video);
        img = mask.*double(img(:,:,1));   %just take first color matrix
        
        %make sure our fly will be a bright spot
        if (strcmp(opt1,'neg') == 1)   
            img = -img+flipper;
        elseif (strcmp(opt1,'pos') == 1)
            img = img;
        end
         
        % Calculate square of difference image (to reduce noise) and track
        temp = (img-background).^2;
              
        [x_raw y_raw bodyline sqr]=FindFly(temp,sqrsize); 
        
        
        %Just in case of a false track due to objects in the video that are
        %brighter than fly and vary over time.  We assume that these are
        %very far away from the fly.  One may want to change the threshold
        %from 100 pixels to something smaller, if need be.  Any reasonably
        %large number should do, and if you choose an incorrect threshold,
        %you can use the CleanData function to fix this.
        if (frame~=StartFrame)
            if (sqrt((x_raw-x(end))^2+(y_raw-y(end))^2) >= 100)
                x_raw = x(end);
                y_raw = y(end);
            end
        end
        
        x = [x x_raw];
        y = [y y_raw];
        orientUHP = [orientUHP bodyline(1)];
        orientLHP = [orientLHP bodyline(2)];
% %         
%         figure(2)
%         imagesc(img)
%         colormap gray
%         hold on
%         plot(x_raw,y_raw, '.r')
% %        
%         
        %grab whatever the background was before the fly got there.
        temp = background(sqr(1):sqr(2),sqr(3):sqr(4));
        %update
        background = alpha*background+(1-alpha)*img;
        %make sure we only update image outside of fly
        background(sqr(1):sqr(2),sqr(3):sqr(4)) = temp;
       
        worked = next(video);
        if (worked == 0)
            break
        end
end
disp('Finished tracking.')
disp('Have a lovely day.')
toc

%Now we have to make sure that the y vector is based on a coordinate
%system with the origin in the lower left hand corner of the square, 
%for plotting purposes.
y = height-y;

orientation = [orientUHP;orientLHP];

return; 
    