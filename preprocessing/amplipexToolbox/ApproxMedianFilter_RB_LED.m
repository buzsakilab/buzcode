function whl = ApproxMedianFilter_RB_LED(file)

% This code tests an approximate median filter for separating out the red
% and blue leds mounted on top of the subjects head
% I was using the DataMax system, which was synced using a red LED. Hence,
% I am tracking the x,y position of the animal as well as when the sync
% light comes on.


%file       = 'Mouse12-120808-01.mpg';
%[~,name,~] = fileparts(file);
% Create Text file for writing out LED locations and sync trigger
%fid        = fopen(sprintf('%s.whl',name),'w+');
% Creater readerobj of file
readerobj  = VideoReader(file);
width      = readerobj.Width;
height     = readerobj.Height;
threshF    = 110;    % threshold on foreground pixel intensity
% Initial frame
Fint       = 1;

% Initialize grid for locating centroid of LED
[X,Y] = meshgrid(1:width,1:height);

% Initialize background as a grayscale image of the first frame
bg_bw     = rgb2gray(read(readerobj,Fint));

% Initialize foreground image and difference image
fg          = zeros(size(bg_bw));
fr_diff     = zeros(size(bg_bw));

% Initialize color mask
mask  = zeros(height,width,3,'uint8');

% Initialize fr in case first frame reading returns error
fr  = zeros(height,width,3,'uint8');

% Initialize whl matrix
whl = zeros(readerobj.NumberOfFrames,4);

for i = Fint:readerobj.NumberOfFrames
    %fprintf('%i',i);
    
    % Access frame of interest - if error (mostly for the last frames,
    % don't know why), will analyze previous frame...
    try
        fr    = read(readerobj,i);
    end
    % convert frame to grayscale
    fr_bw = rgb2gray(fr);
%     keyboard
    %%% Label Color Mask
    label         = repmat(logical(fr_bw>threshF),[1 1 3]);
    mask(label)   = fr(label);
    mask(~label)  = 0;
    
    %%% Find centroid of remaining pixels %%%
    bw_mask = rgb2gray(mask);
    [L,Nl]  = bwlabel(double(bw_mask));
    ids     = [];
    
    % Initialize red and blue LEDS as missing (=> -1)
    Rr = [-1 -1];
    Br = [-1 -1];
    
    % Update centroid if connected components are found
    if Nl > 0
        a = regionprops(L,'PixelList');
        
        for ii = 1:Nl
            % Skip single pixels
            if size(a(ii).PixelList,1)<2
                continue
            end
            
            % Access color information and calculate means
            R  = mask(a(ii).PixelList(:,2),a(ii).PixelList(:,1),1);
            mR = mean(R(R>0));
            B  = mask(a(ii).PixelList(:,2),a(ii).PixelList(:,1),3);
            mB = mean(B(B>0));
            
            if mR > mB && mR > 200
                ids = [ids; a(ii).PixelList];
                Rr   = round(mean(a(ii).PixelList,1));
            elseif mR < mB && mB > 200
                ids = [ids; a(ii).PixelList];
                Br   = round(mean(a(ii).PixelList,1));
            else
                continue
            end
            
        end
    end
    
    whl(i,:) = [Rr(1),Rr(2),Br(1),Br(2)];
  
    % End processing time
    if mod(i,100)==0 & i<1000
        h = waitbar(i/readerobj.NumberOfFrames);
        if 1
            ixr = whl(i-99:i,1);
            ixr = ixr>-1;
            ixb = whl(i-99:i,3);
            ixb = ixb>-1;
            figure(1),clf,
            subplot(3,1,1);imshow(uint8(fr));title('If looks wrong, Ctrl+C and change threshold \!!')
            subplot(3,1,2),imshow(uint8(bw_mask));title('Mask')
            subplot(3,1,3)
                plot(whl(ixr,1),whl(ixr,2),'r')
                hold on
                plot(whl(ixb,3),whl(ixb,4),'b')
            fprintf('Detection of Red LED failed %i times (%i times for the blue LED) in the last 100 frames\n',sum(~ixr),sum(~ixb));   
        end
    elseif mod(i,1000)==0
        h = waitbar(i/readerobj.NumberOfFrames);
    end
  
end
close(h)