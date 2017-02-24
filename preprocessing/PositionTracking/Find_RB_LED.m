function whl = Find_RB_LED(file)

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
threshF    = 150;    % threshold on foreground pixel intensity
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
rmask  = zeros(height,width,3,'uint8');
bmask  = zeros(height,width,3,'uint8');

% Initialize fr in case first frame reading returns error
fr  = zeros(height,width,3,'uint8');

% Initialize whl matrix
whl = zeros(readerobj.NumberOfFrames,4);

for i = Fint:readerobj.NumberOfFrames
%     % Initialize red and blue LEDS as missing (=> -1)
    Red = [-1 -1];
    Blue = [-1 -1];
    
    %fprintf('%i',i);
    
    % Access frame of interest - if error (mostly for the last frames,
    % don't know why), will analyze previous frame...
    try
        fr    = read(readerobj,i);
    end
    
    % Spatial filter
    filt = fspecial('average',50);
    frs = imfilter(fr,filt);
    fr2 = fr-frs;

    filt = fspecial('average',10);
    fr3 = imfilter(fr2, filt);

    % convert frame to thresholdable
%     fr_bw = rgb2gray(fr);
    fr_r = double(fr2(:,:,1));
    fr_b = double(fr2(:,:,3));
    template = zthreshimage(rgb2gray(fr2),0.1);%pixels that aren't too dim after filtering... ie worth paying attention to
    fr_r_ok = fr_r;
    fr_b_ok = fr_b;
    fr_b_ok(find(~template)) = 0;
    fr_r_ok(find(~template)) = 0;
    
    rbdiff = fr_r-fr_b;
    brdiff = fr_b-fr_r;
    %would have to brightness threshold first
    rbratio = log(fr_r_ok./fr_b_ok);
    brratio = log(fr_b_ok./fr_r_ok);
    
    % bright enough in each color
    zbright_r = zthreshimage(fr_r,1);
    zbright_b = zthreshimage(fr_b,1);
    absbright_r = threshimage(fr_r,60);
    absbright_b = threshimage(fr_b,60);
    
    % color selective enough by diff
    zdiffsel_r = zthreshimage(rbdiff,3);% diff above z score
    zdiffsel_b = zthreshimage(brdiff,3);
    absdiffsel_r = threshimage(rbdiff,20);% diff above absolute level 
    absdiffsel_b = threshimage(brdiff,20);

    % color selective enough by ratio
    zratsel_r = zthreshimage(rbratio,1);% diff above z score
    zratsel_b = zthreshimage(brratio,1);
    
%     ratsel_r = im2bw(rbratio,.51);% above 55st %ile ratio
%     ratsel_b = im2bw(brratio,.51);% above 51st %ile ratio
    
    % bright and selective enough
%     combo_r = diffsel_r.*ratsel_r.*bright_r;
%     combo_b = diffsel_b.*ratsel_r.*bright_b;
    combo_r = zbright_r.*absbright_r.*zdiffsel_r.*absdiffsel_r.*zratsel_r;
    combo_b = zbright_b.*absbright_b.*zdiffsel_b.*absdiffsel_b.*zratsel_b;
    
% keep only the largest object in each thresholded image
    % blue
    [l,numobjs] = bwlabel(combo_b);
    if numobjs>0
        numeach = zeros(1,numobjs);
        for a = 1:numobjs
            numeach(a) = sum(l(:)==a);
        end
        [~,bgroup] = max(numeach);
        pl = regionprops(l,'PixelList');
        Blue = round(mean(pl(bgroup).PixelList,1));
    end
    % red
    [l,numobjs] = bwlabel(combo_r);
    if numobjs>0
        numeach = zeros(1,numobjs);
        for a = 1:numobjs
            numeach(a) = sum(l(:)==a);
        end
        [~,rgroup] = max(numeach);
        pl = regionprops(l,'PixelList');
        Red = round(mean(pl(rgroup).PixelList,1));
    end
    
    whl(i,:) = [Red(1),Red(2),Blue(1),Blue(2)];
  
    % End processing time, now outputting/plotting
    if mod(i,100)==0 & i<1000
        ixr = whl(i-99:i,1);
        ixr = ixr>-1;
        ixb = whl(i-99:i,3);
        ixb = ixb>-1;
        ok = ixr.*ixb;
        figure(1),clf,
        subplot(3,1,1);imshow(uint8(fr));title('If looks wrong, Ctrl+C and change threshold \!!')

        subplot(3,3,5);imagesc(cat(3,combo_r,zeros(size(combo_r)),combo_b));
        xlim([0 size(fr,1)])
        ylim([0 size(fr,2)])
        axis tight

        subplot(3,3,8)
            plot(whl(i-99:i,1),whl(i-99:i,2),'r.')
            hold on
            plot(whl(i-99:i,3),whl(i-99:i,4),'b.')
            mx = mean([whl(i-100+find(ok),1) whl(i-100+find(ok),3)],2);                
            my = mean([whl(i-100+find(ok),2) whl(i-100+find(ok),4)],2);                
            plot(mx,my,'Marker','.','color',[0.5 0.5 0.5],'LineStyle','none')
            xlim([0 size(fr,1)])
            ylim([0 size(fr,2)])
            axis ij
        fprintf('Detection of Red LED failed %i times (%i times for the blue LED) in the last 100 frames\n',sum(~ixr),sum(~ixb));   
    end
    if i == 1;
        h = waitbar(i/readerobj.NumberOfFrames);
    else
        waitbar(i/readerobj.NumberOfFrames,h);
    end
  
end
close(h)